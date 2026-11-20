# Two-group randomization (permutation) test with effect sizes and intervals.
#' Randomization (Permutation) Test for Two Independent Groups
#'
#' Compares two independent groups by referring the observed statistic to
#' the distribution of that same statistic over reassignments of the
#' observed scores to the two groups. That reference distribution, not a
#' normal or a \emph{t} distribution, supplies the \emph{p}-value, so the
#' test needs no assumption about the shape of the population. Alongside
#' the test, the function reports the effect sizes that answer the
#' question the test only screens: the mean difference and its
#' randomization-based interval, the standardized mean difference with a
#' noncentral \emph{t} interval, the common language effect size, and
#' Cliff's delta.
#'
#' @param x A formula of the form \code{y ~ group}, or a numeric response
#'   vector to be paired with \code{group}. Leave \code{NULL} when the two
#'   samples are supplied through \code{group_1} and \code{group_2}.
#' @param group A grouping variable the same length as \code{x}, with
#'   exactly two levels after unused levels are dropped. Ignored when
#'   \code{x} is a formula.
#' @param data An optional \code{data.frame} in which to find the variables
#'   named in the formula.
#' @param group_1,group_2 The two samples supplied directly as numeric
#'   vectors, an alternative to the formula and response-plus-grouping
#'   interfaces. Lengths need not be equal.
#' @param statistic One of \code{"mean"} (default) or \code{"t"}.
#'   \code{"mean"} uses the difference in group means. \code{"t"} uses the
#'   studentized (Welch) statistic, which divides that difference by its
#'   separate-variances standard error and is the better choice when the
#'   groups differ in size and may differ in variance (see Details). The
#'   statistic is refused when both observed groups are constant at
#'   different values (exactly, or to within the tolerance of the
#'   comparison), where its standard error is 0 (see \emph{When the
#'   studentized statistic is undefined}); a reassignment with that
#'   property in the reference set of other data is ordered as \code{Inf}
#'   or \code{-Inf}, and when every score is equal the statistic is 0.
#' @param alternative One of \code{"two_sided"} (default; the base-R
#'   spelling \code{"two.sided"} is accepted as an alias), \code{"less"},
#'   or \code{"greater"}. The direction refers to the first group minus the
#'   second, the same orientation \code{\link[stats]{t.test}} uses.
#' @param exact Logical. If \code{NULL} (default), the test enumerates
#'   every reassignment when \code{choose(N, n_1)} is at most 50,000 and
#'   samples reassignments otherwise. If \code{TRUE}, enumeration is forced
#'   (refused above 1,000,000 reassignments). If \code{FALSE}, Monte Carlo
#'   is forced.
#' @param n_resamples A single whole number of at least 1: the number of
#'   randomly drawn reassignments when enumeration is not used. Default
#'   \code{10000L}.
#' @param seed Optional integer seed for the Monte Carlo branch. Default
#'   \code{NULL}, under which the reassignments are drawn from, and
#'   advance, the current random number stream; supply an integer for
#'   reproducibility. When a seed is supplied the RNG state in place
#'   before the call is restored on exit.
#' @param conf_level Confidence level for every interval reported, the
#'   inverted randomization interval included. Default \code{0.95}.
#' @param shift_ci Logical. Compute the randomization-based interval for
#'   the shift by inverting the test? Default \code{TRUE}. Setting it to
#'   \code{FALSE} reports the two endpoints as \code{NA} and skips the
#'   inversion, which is the expensive part of the call.
#'
#' @return A \code{data.frame} with a \code{term} column and a numeric
#'   \code{value} column, in three blocks.
#'
#'   The test: \code{mean_difference} (first group minus second),
#'   \code{statistic} (the statistic actually referred to the reference
#'   distribution, always finite), \code{p_value}, and \code{p_value_se}
#'   (the Monte Carlo standard error of the \emph{p}-value, \code{NA}
#'   under exact enumeration, which has no Monte Carlo error).
#'
#'   The intervals and effect sizes: \code{shift_ci_lower} and
#'   \code{shift_ci_upper} (the randomization interval for the shift,
#'   obtained by inverting the test); \code{normal_theory_ci_lower} and
#'   \code{normal_theory_ci_upper} (Welch's \emph{t} interval on the same
#'   mean difference, reported for contrast); \code{smd} with
#'   \code{smd_ci_lower} and \code{smd_ci_upper} from
#'   \code{\link{ci_smd}}; \code{cles} with \code{cles_ci_lower} and
#'   \code{cles_ci_upper} from \code{\link{cles}}; and
#'   \code{cliff_delta} with \code{cliff_delta_ci_lower} and
#'   \code{cliff_delta_ci_upper} from \code{\link{cliff_delta}}. When
#'   the pooled within-group variance is zero (each group constant) the
#'   standardized mean difference is undefined: the six rows \code{smd},
#'   \code{smd_ci_lower}, \code{smd_ci_upper}, \code{cles},
#'   \code{cles_ci_lower}, and \code{cles_ci_upper} are \code{NA},
#'   a warning of class \code{dmar_smd_undefined} says so, and the test,
#'   the shift interval, and Cliff's delta are reported as usual.
#'
#'   The design: \code{n_1}, \code{n_2}, \code{N}, \code{n_evaluated} (how
#'   many reassignments were actually used), and \code{exact} (1 if every
#'   reassignment was enumerated, 0 if they were sampled).
#'
#'   Non-numeric information travels on attributes rather than in the
#'   \code{value} column: \code{statistic_name}, \code{method}
#'   (\dQuote{exact enumeration} or \dQuote{Monte Carlo}),
#'   \code{alternative}, \code{group_labels}, \code{response_name},
#'   \code{group_name}, \code{seed}, \code{observed_statistic}, and
#'   \code{reference_distribution}, the vector of statistics over the
#'   reassignments that \code{\link{plot_randomization_test}} draws, which
#'   can hold \code{Inf} and \code{-Inf} under \code{statistic = "t"}
#'   (see \code{statistic}). \code{reference_extreme} is a logical vector
#'   of the same length that marks the reassignments the \emph{p}-value
#'   counted as at least as extreme as the observed one.
#'
#'   No table is returned when \code{statistic = "t"} is refused; the
#'   condition, of class \code{dmar_studentized_statistic_undefined},
#'   carries \code{reason} (\code{"constant_groups"}),
#'   \code{alternative}, \code{n_1}, \code{n_2}, and
#'   \code{p_value_mean}, the exact \emph{p}-value
#'   \code{statistic = "mean"} gives for the same data and alternative.
#'
#' @details
#' \strong{What the randomization distribution is.} Suppose \emph{N}
#' participants were randomly assigned, \eqn{n_1} to one condition and
#' \eqn{n_2} to the other. Under the null hypothesis that the condition a
#' participant received made no difference to that participant's score,
#' each score would have been the same number no matter which group the
#' participant landed in. The assignment actually used was one draw from
#' the \eqn{\binom{N}{n_1}} assignments the randomization could equally
#' well have produced, so every one of those assignments was equally
#' likely, and each of them yields a value of the test statistic. Those
#' values are the randomization distribution. The \emph{p}-value is the
#' proportion of them at least as extreme as the value the experiment
#' actually produced.
#'
#' \strong{Why there is no normality assumption.} Nothing in that argument
#' mentions a population, a normal curve, or a sampling model. The
#' probability comes from the coin flips the experimenter performed, which
#' are known exactly because the experimenter performed them. This is the
#' inferential logic Fisher (1935) used to introduce experimental design,
#' and it is where Chapter 1 of Maxwell, Delaney, and Kelley (2027) starts,
#' for the same reason: the validity of the test rests on the randomization
#' rather than on assumptions a data analyst cannot check.
#'
#' \strong{What the test does and does not license.} A small
#' \emph{p}-value says the observed separation between the groups would
#' rarely arise from reassignment alone, which is evidence that the
#' assignment mattered. It does not say how much it mattered, and with a
#' large \emph{N} an uninteresting difference will produce a small
#' \emph{p}-value. It also does not, by itself, license generalization
#' beyond the participants at hand: randomization licenses a causal claim
#' about these units, while generalization to a population is a separate
#' argument that rests on how the units were recruited. That is why this
#' function reports effect sizes with intervals rather than a
#' \emph{p}-value alone.
#'
#' \strong{Why the studentized statistic.} With \eqn{n_1 = n_2} the
#' separate-variances standard error equals the pooled one, and the pooled
#' \emph{t} is an increasing function of the mean difference over the
#' reassignments, so the two statistics give the same \emph{p}-value in
#' every balanced design. They differ only when the groups are unequal in
#' size. Reassigning scores between groups of unequal size
#' mixes the two variances in proportions that the observed assignment does
#' not have, so the reference distribution for the raw mean difference is
#' built under a null that is false in a second way, and the test's actual
#' Type I error rate drifts away from the nominal level. The studentized
#' statistic rescales each reassignment by its own separate-variances
#' standard error, which removes most of that drift and remains
#' asymptotically valid under heteroscedasticity (Janssen, 1997; Neuhaus,
#' 1993). Use \code{statistic = "t"} when the group sizes differ and
#' unequal variances are plausible.
#'
#' \strong{When the studentized statistic is undefined.} The studentized
#' statistic is a ratio, the difference in means over its standard error,
#' and it has a value only when that standard error is positive. If each
#' group is constant and the two constants differ (exactly, or to within
#' the tolerance of the comparison), the standard error of the observed
#' data is 0 and \code{statistic = "t"} is refused with a condition of
#' class \code{dmar_studentized_statistic_undefined}. The observed data
#' are then at the boundary of the reference set: among the reassignments
#' only the observed one and, with equal group sizes, the one that
#' exchanges the groups have a zero standard error. The refusal concerns
#' the statistic and not the evidence. \code{statistic = "mean"} is
#' defined for every arrangement and gives the smallest \emph{p}-value the
#' design allows in the direction of the data,
#' \eqn{1/\binom{N}{n_1}} one-sided, and two-sided
#' \eqn{2/\binom{N}{n_1}} with equal group sizes and
#' \eqn{1/\binom{N}{n_1}} otherwise; the message states the value, and
#' the condition carries it as \code{p_value_mean}. The refusal is made
#' before any reassignment is enumerated or drawn, so a refused Monte
#' Carlo call leaves the random number stream where it was. With
#' \code{statistic = "mean"} the same data return the table, with the
#' warning of class \code{dmar_smd_undefined} described under Value.
#'
#' When the observed statistic is defined, some arrangements in the
#' reference set can still have a zero standard error: with a two-valued
#' outcome, the reassignment that puts every high score in one group.
#' Those arrangements are kept and ordered as the limit of the statistic,
#' \code{Inf} or \code{-Inf} with the sign of their mean difference, so
#' they count as the most extreme arrangements, which is where the mean
#' difference also places them. Refusing them would refuse ordinary binary
#' and coarsely measured data. For \code{group_1 = c(1, 1, 1, 0)} against
#' \code{group_2 = c(0, 0, 0, 1)} the reference distribution of the
#' studentized statistic holds one \code{-Inf} and one \code{Inf} among
#' its 70 values, and the \emph{p}-value, 34/70, is the one the difference
#' in means gives. When every score is equal, every arrangement is the
#' same data set; the statistic is reported as 0 and the \emph{p}-value
#' is 1.
#'
#' \strong{Ties.} A reassignment whose statistic equals the observed one
#' counts as at least as extreme. Two statistics that are equal in exact
#' arithmetic can differ in their last binary digits as computed, so the
#' comparison allows for the rounding of the arithmetic, with a bound
#' proportional to the size of the scores. Where the bound cannot separate
#' two statistics they are counted as tied, so the \emph{p}-value can be
#' above the exact one and is not below it. For the difference in means
#' that happens only when two mean differences agree to within about 1e-14
#' of the largest absolute score, as when scores agree to about fourteen
#' significant digits or one score is about 1e14 times the rest. For the
#' studentized statistic it can also happen when one score is about ten
#' million times the spread of the rest, or when every score lies within
#' about a millionth of one of two values.
#'
#' \strong{Exact or Monte Carlo.} When \code{choose(N, n_1)} is at most
#' 50,000 every reassignment is enumerated and the \emph{p}-value is exact:
#' it is a count divided by a known total, with no approximation anywhere.
#' Above that threshold \code{n_resamples} reassignments are drawn at
#' random and the \emph{p}-value is
#' \eqn{(r + 1) / (m + 1)}, where \emph{r} counts the sampled reassignments
#' at least as extreme as the observed one and \emph{m} is
#' \code{n_resamples}. Adding one to each part counts the observed
#' assignment, which is itself a legitimate reassignment; without it a
#' \emph{p}-value of exactly zero could be reported for a hypothesis the
#' data cannot rule out, and the test would be anticonservative (Phipson &
#' Smyth, 2010). The reported \code{p_value_se} is
#' \eqn{\sqrt{\hat p (1 - \hat p) / m}}, the standard error of the
#' resampling itself. It describes how much the \emph{p}-value would move
#' if the reassignments were drawn again, not how much it would move in a
#' new experiment. Raising \code{n_resamples} shrinks it at the usual
#' \eqn{1/\sqrt{m}} rate.
#'
#' \strong{The randomization interval, and how it differs from the normal
#' theory one.} Suppose the treatment adds a constant \eqn{\delta} to every
#' score it touches. Subtracting \eqn{\delta} from each first-group score
#' should then leave scores that are exchangeable across groups, so the
#' randomization test applied to the subtracted data is a test of
#' \eqn{H_0\!: \mathrm{shift} = \delta}. The set of \eqn{\delta} for which
#' that test does not reject at level \eqn{1 - } \code{conf_level} is a
#' confidence interval for the shift, and it is reported as
#' \code{shift_ci_lower} and \code{shift_ci_upper}. Inverting a test
#' this way is the general recipe (Ernst, 2004); the endpoints are located
#' by bisection on the \emph{p}-value, using the same reassignments
#' throughout so the interval and the test agree.
#'
#' The contrast with \code{normal_theory_ci_lower} and
#' \code{normal_theory_ci_upper}, which are Welch's \emph{t} limits on
#' the same mean difference, is worth reading whenever both are printed.
#' The randomization interval is exactly the set of shifts the test being
#' run does not reject, so the test and the interval can never disagree.
#' The normal theory interval instead assumes the sampling distribution of
#' the mean difference has a known shape; it is smooth, symmetric about the
#' point estimate, and can extend past the range the data can support.
#' The randomization interval is discrete, need not be symmetric, and in a
#' very small design is unbounded, a correct statement of how little
#' information the design carries rather than a defect:
#' with three observations per group the smallest attainable two-sided
#' \emph{p}-value is 2/20 = 0.10, so no shift can be rejected at the 5\%
#' level and the 95\% interval is the whole real line. The randomization
#' interval also inherits the shift model, so it answers a narrower
#' question than the test does: the test needs only exchangeability, while
#' the interval needs the treatment to move every score by the same amount.
#'
#' \strong{Effect sizes.} Every effect size reported here comes from the
#' package function that owns it, so the numbers match a direct call. The
#' one exception is a pooled within-group variance of zero, where the
#' standardized mean difference and the common language effect size are
#' reported as \code{NA} (see Value).
#' \code{\link{smd}} and \code{\link{ci_smd}} supply the standardized mean
#' difference and its noncentral \emph{t} interval; \code{\link{cles}}
#' supplies the common language effect size, the probability that a
#' randomly drawn score from the first group exceeds one from the second,
#' by transforming those limits through \eqn{\Phi(\cdot/\sqrt 2)}; and
#' \code{\link{cliff_delta}} supplies Cliff's delta with its consistent
#' interval. The standardized mean difference and the common language
#' effect size are normal theory quantities, so their intervals lean on the
#' assumption the test itself avoids. Cliff's delta does not: it is a
#' function of the ordering of the observations alone, which makes it the
#' natural effect size companion to a randomization test. Reporting all
#' three lets a reader see whether the distribution-free and normal theory
#' summaries tell the same story.
#'
#' \strong{Ranks give the Wilcoxon test.} Replacing the scores by their
#' ranks and running this test with \code{statistic = "mean"} reproduces
#' the exact Wilcoxon rank sum test, since the rank sum is a monotone
#' function of the difference in mean ranks. That equivalence is a useful
#' check and a reminder of what the rank test is: a randomization test on
#' transformed data.
#'
#' @references
#' Edgington, E. S., & Onghena, P. (2007). \emph{Randomization tests}
#'   (4th ed.). Chapman & Hall/CRC.
#'
#' Ernst, M. D. (2004). Permutation methods: A basis for exact inference.
#'   \emph{Statistical Science, 19}(4), 676--685.
#'   \doi{10.1214/088342304000000396}
#'
#' Fisher, R. A. (1935). \emph{The design of experiments}. Oliver & Boyd.
#'
#' Janssen, A. (1997). Studentized permutation tests for non-i.i.d.
#'   hypotheses and the generalized Behrens-Fisher problem.
#'   \emph{Statistics & Probability Letters, 36}(1), 9--21.
#'   \doi{10.1016/S0167-7152(97)00043-6}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 1 on the logic of randomization
#'   and the randomization test.)
#'
#' Neuhaus, G. (1993). Conditional rank tests for the two-sample problem
#'   under random censorship. \emph{The Annals of Statistics, 21}(4),
#'   1760--1779. \doi{10.1214/aos/1176349396}
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation \emph{p}-values should
#'   never be zero: Calculating exact \emph{p}-values when permutations are
#'   randomly drawn. \emph{Statistical Applications in Genetics and
#'   Molecular Biology, 9}(1), Article 39. \doi{10.2202/1544-6115.1585}
#'
#' Pitman, E. J. G. (1937). Significance tests which may be applied to
#'   samples from any populations. \emph{Supplement to the Journal of the
#'   Royal Statistical Society, 4}(1), 119--130.
#'
#' @seealso \code{\link{plot_randomization_test}} for the figure that shows
#'   the reference distribution, the observed statistic, and the rejection
#'   region; \code{\link{randomization_test_paired}} for the sign-flip
#'   sibling used with paired observations;
#'   \code{\link[stats]{t.test}} and \code{\link[stats]{wilcox.test}} for
#'   the parametric and rank-based alternatives; \code{\link{smd}},
#'   \code{\link{ci_smd}}, \code{\link{cles}}, and
#'   \code{\link{cliff_delta}} for the effect sizes reported here.
#'
#' @examples
#' # 1. Ten observations, so every one of the choose(10, 5) = 252
#' #    reassignments is enumerated and the p-value is exact.
#' treatment <- c(80, 84, 79, 88, 83)
#' control   <- c(72, 75, 68, 81, 74)
#' randomization_test(group_1 = treatment, group_2 = control)
#'
#' # 2. The studentized statistic. The groups are the same size here, so it
#' #    gives the same p-value as the difference in means; the two differ
#' #    when the groups are unequal in size.
#' randomization_test(group_1 = treatment, group_2 = control,
#'                    statistic = "t")
#'
#' # 3. Formula interface: weekly drinking in the two comparable arms of
#' #    the drinks_trial data, a right-skewed outcome, which is exactly
#' #    where a distribution-free test earns its keep. With 37 and 32
#' #    participants there are far too many reassignments to enumerate, so
#' #    they are drawn at random and the p-value carries a Monte Carlo
#' #    standard error. n_resamples = 2000 keeps the example quick; a
#' #    reported p-value deserves the default of 10,000.
#' cra <- droplevels(subset(drinks_trial, treatment != "CRA + Disulfiram"))
#' set.seed(113)
#' randomization_test(drinks_per_week ~ treatment, data = cra,
#'                    n_resamples = 2000, seed = 113)
#'
#' # 4. A one-sided test, and the one-sided interval that goes with it.
#' randomization_test(group_1 = treatment, group_2 = control,
#'                    alternative = "greater")
#'
#' # 5. On ranks, the test is the exact Wilcoxon rank sum test.
#' y <- c(treatment, control)
#' g <- rep(c("treatment", "control"), each = 5)
#' res <- randomization_test(rank(y), g)
#' res$value[res$term == "p_value"]
#' wilcox.test(treatment, control, exact = TRUE)$p.value
#'
#' # 6. Both observed groups are constant, at different values. The
#' #    difference in means is defined, and gives the smallest p-value
#' #    groups of 4 and 4 allow. The warnings say that the standardized
#' #    mean difference is undefined and that the groups do not overlap.
#' #    With statistic = "t" the same call is refused, with a message that
#' #    says why: the standard error is 0, so the studentized statistic is
#' #    undefined.
#' randomization_test(group_1 = rep(5, 4), group_2 = rep(3, 4))
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest nonparametric
#'
#' @family hypothesis tests
#'
#' @export

randomization_test <- function(x = NULL, group = NULL, data = NULL,
                               group_1 = NULL, group_2 = NULL,
                               statistic = c("mean", "t"),
                               alternative = c("two_sided", "less", "greater"),
                               exact = NULL,
                               n_resamples = 10000L,
                               seed = NULL,
                               conf_level = 0.95,
                               shift_ci = TRUE) {
  .dmar_check_conf_level(conf_level)
  statistic   <- match.arg(statistic)
  alternative <- .match_alternative(alternative)

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1)
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  n_resamples <- .dmar_check_count(n_resamples, "n_resamples")

  # ---------- Resolve the three input interfaces ----------
  # Two samples given directly, a formula, or a response with a grouping
  # variable. Two numeric vectors of the same length are ambiguous between
  # "two samples" and "response plus a numeric grouping code", so the two
  # sample form is reached only through the named group_1 / group_2
  # arguments and never by guessing.
  response_name <- NULL
  group_name    <- NULL

  if (!is.null(group_1) || !is.null(group_2)) {
    if (!is.null(x) || !is.null(group))
      .dmar_refuse_argument(
        "Supply the data either as 'group_1' and 'group_2' or as a formula (or a response with 'group'), not both.",
        arg = c("group_1", "group_2", "group"))
    if (is.null(group_1) || is.null(group_2))
      .dmar_refuse_argument(
        "Both 'group_1' and 'group_2' are needed when the samples are supplied directly.",
        arg = c("group_1", "group_2"))
    if (!is.numeric(group_1) || !is.numeric(group_2))
      .dmar_refuse_argument(
        "'group_1' and 'group_2' must both be numeric vectors.",
        arg = c("group_1", "group_2"))
    group_1 <- group_1[!is.na(group_1)]
    group_2 <- group_2[!is.na(group_2)]
    y <- c(group_1, group_2)
    g <- factor(rep(c("group_1", "group_2"),
                    c(length(group_1), length(group_2))),
                levels = c("group_1", "group_2"))
  } else if (inherits(x, "formula")) {
    mf <- if (is.null(data)) stats::model.frame(x) else
      stats::model.frame(x, data = data)
    if (ncol(mf) != 2L)
      .dmar_refuse_argument(
        "The formula must have the form 'y ~ group': one response and one grouping variable.",
        arg = NULL)
    y <- mf[[1L]]
    g <- mf[[2L]]
    response_name <- names(mf)[1L]
    group_name    <- names(mf)[2L]
  } else if (is.numeric(x) && !is.null(group)) {
    if (length(x) != length(group))
      .dmar_refuse_argument(
        "'x' and 'group' must be the same length. To supply two samples directly, use 'group_1' and 'group_2'.",
        arg = c("x", "group", "group_1"))
    y <- x
    g <- group
  } else {
    .dmar_refuse_argument(
      "Supply a formula 'y ~ group' (with 'data'), a numeric response 'x' with a 'group' variable, or the two samples as 'group_1' and 'group_2'.",
      arg = c("data", "x", "group"))
  }

  if (!is.numeric(y)) .dmar_refuse_argument(
    "The response must be numeric.", arg = NULL)
  if (!is.factor(g)) g <- factor(g)
  ok <- !is.na(y) & !is.na(g)
  y  <- y[ok]
  g  <- droplevels(g[ok])
  if (!all(is.finite(y)))
    .dmar_refuse_argument("The response must be finite.", arg = NULL)

  lev <- levels(g)
  if (length(lev) != 2L)
    .dmar_refuse_argument(c(
      "randomization_test() compares exactly two groups; the grouping variable has ",
      length(lev), " level(s) with data."), arg = NULL)

  in_1 <- g == lev[1L]
  n_1  <- sum(in_1)
  n_2  <- sum(!in_1)
  if (n_1 < 2L || n_2 < 2L)
    .dmar_refuse_argument(
      "Each group needs at least 2 observations.", arg = NULL)
  N <- n_1 + n_2

  y_1 <- y[in_1]
  y_2 <- y[!in_1]
  mean_difference <- mean(y_1) - mean(y_2)
  v_1 <- stats::var(y_1)
  v_2 <- stats::var(y_2)
  se_welch <- sqrt(v_1 / n_1 + v_2 / n_2)
  df_welch <- if (se_welch > 0)
    (v_1 / n_1 + v_2 / n_2)^2 /
      ((v_1 / n_1)^2 / (n_1 - 1) + (v_2 / n_2)^2 / (n_2 - 1)) else NA_real_

  # ---------- Enumerate or sample the reassignments ----------
  n_possible <- choose(N, n_1)
  use_exact  <- if (is.null(exact)) n_possible <= .rt_exact_threshold else isTRUE(exact)
  if (use_exact && n_possible > .rt_exact_cap)
    .dmar_refuse_argument(c(
      "Exact enumeration would require ", format(n_possible, big.mark = ",", scientific = FALSE),
      " reassignments, above the ", format(.rt_exact_cap, big.mark = ",", scientific = FALSE),
      " cap. Use exact = FALSE for the Monte Carlo test."), arg = NULL)

  obs_idx_1 <- which(in_1)
  mem       <- as.numeric(seq_len(N) %in% obs_idx_1)

  # The observed assignment goes through the same arithmetic as the
  # reassignments, and it goes first: when both observed groups are
  # constant at different values the studentized statistic has a zero
  # standard error and no value, and the call is refused before a
  # reassignment is enumerated or drawn. The observed standard error does
  # not depend on the shift, so a call that passes has a defined observed
  # statistic at every shift of the interval search.
  w_obs <- .rt_ind_workspace(y, mem, matrix(obs_idx_1, ncol = 1L), statistic)
  obs_0 <- .rt_ind_null(w_obs, 0)
  if (statistic == "t" && is.infinite(obs_0$statistic))
    .rt_refuse_constant_groups(y_1, y_2, alternative, n_possible, use_exact,
                               forced = identical(exact, FALSE) &&
                                 n_possible <= .rt_exact_threshold,
                               n_resamples = n_resamples)

  if (use_exact) {
    perm   <- utils::combn(N, n_1)
    n_eval <- ncol(perm)
  } else {
    .dmar_local_seed(seed)
    perm   <- matrix(replicate(n_resamples, sample.int(N, n_1)), nrow = n_1)
    n_eval <- n_resamples
  }

  # Each statistic comes with its rounding bound, and a reference value
  # within the sum of its own bound and the observed value's bound of the
  # observed statistic is a tie.
  w_null <- .rt_ind_workspace(y, mem, perm, statistic)
  null_0 <- .rt_ind_null(w_null, 0)
  T_null <- null_0$statistic
  T_obs  <- obs_0$statistic

  flags <- .rt_is_extreme(T_null, T_obs, alternative,
                          null_0$tol + obs_0$tol)
  n_extreme <- sum(flags)
  if (use_exact) {
    p_value <- n_extreme / n_eval
    p_se    <- NA_real_
  } else {
    p_value <- (1 + n_extreme) / (1 + n_eval)
    p_se    <- sqrt(p_value * (1 - p_value) / n_eval)
  }

  # ---------- Invert the test for the shift interval ----------
  alpha <- 1 - conf_level
  shift_lower <- NA_real_
  shift_upper <- NA_real_
  if (isTRUE(shift_ci)) {
    p_fun <- function(delta) {
      null_d <- .rt_ind_null(w_null, delta)
      obs_d  <- .rt_ind_null(w_obs,  delta)
      r <- .rt_extreme(null_d$statistic, obs_d$statistic, alternative,
                       null_d$tol + obs_d$tol)
      if (use_exact) r / n_eval else (1 + r) / (1 + n_eval)
    }
    # A step of one standard error brackets quickly; the tolerance is small
    # relative to that step, and the bisection halts on it.
    step0 <- max(se_welch, sqrt(.Machine$double.eps) * max(1, abs(mean_difference)))
    tol   <- step0 * 1e-6
    # A shift is rejected when its p-value is at or below alpha. 1 - 0.90 is
    # 0.09999999999999998 in binary arithmetic, below the p-value 1/10 that
    # the test at the 0.10 level rejects, so the comparison is made with
    # alpha raised by a few units of rounding. A p-value is a ratio of
    # counts, and no attainable one lies that close above alpha.
    alpha_shift <- alpha * (1 + 8 * .Machine$double.eps)
    if (alternative != "less")
      shift_lower <- .rt_invert(p_fun, alpha_shift, mean_difference, step0, -1, tol)
    else
      shift_lower <- -Inf
    if (alternative != "greater")
      shift_upper <- .rt_invert(p_fun, alpha_shift, mean_difference, step0, +1, tol)
    else
      shift_upper <- Inf
  }

  # ---------- Normal theory interval on the same mean difference ----------
  if (is.na(df_welch)) {
    nt_lower <- NA_real_
    nt_upper <- NA_real_
  } else if (alternative == "two_sided") {
    half <- stats::qt(1 - alpha / 2, df_welch) * se_welch
    nt_lower <- mean_difference - half
    nt_upper <- mean_difference + half
  } else if (alternative == "less") {
    nt_lower <- -Inf
    nt_upper <- mean_difference + stats::qt(conf_level, df_welch) * se_welch
  } else {
    nt_lower <- mean_difference - stats::qt(conf_level, df_welch) * se_welch
    nt_upper <- Inf
  }

  # ---------- Effect sizes, each from the function that owns it ----------
  # With no within-group variance the standardized mean difference has a
  # zero denominator, and neither it, its noncentral t interval, nor the
  # common language effect size built from them is defined. The test does
  # not depend on any of them, so they are reported as unavailable.
  d <- smd(group_1 = y_1, group_2 = y_2)$value
  if (is.finite(d)) {
    ncp <- d * sqrt(as.numeric(n_1) * n_2 / (n_1 + n_2))
    d_ci <- ci_smd(ncp = ncp, n_1 = n_1, n_2 = n_2, conf_level = conf_level)
    d_lower <- d_ci$value[d_ci$term == "ci_lower"]
    d_upper <- d_ci$value[d_ci$term == "ci_upper"]

    cl <- cles(smd = d, smd_ci_lower = d_lower, smd_ci_upper = d_upper,
               conf_level = conf_level)
    cl_hat   <- cl$value[cl$term == "cles"]
    cl_lower <- cl$value[cl$term == "cles_ci_lower"]
    cl_upper <- cl$value[cl$term == "cles_ci_upper"]
  } else {
    d <- d_lower <- d_upper <- NA_real_
    cl_hat <- cl_lower <- cl_upper <- NA_real_
    .dmar_warn(paste0(
      "The standardized mean difference and the common language effect ",
      "size are undefined because the pooled within-group variance is ",
      "zero; the randomization test and Cliff's delta are reported."),
      "smd_undefined")
  }

  # When the estimated variance of Cliff's delta is 0, cliff_delta() gives
  # a warning of class dmar_interval_degenerate. It is passed on with the
  # same class and elements under this function's name, since this is the
  # function the user called.
  cd <- withCallingHandlers(
    cliff_delta(y_1, y_2, conf_level = conf_level),
    dmar_interval_degenerate = function(w) {
      w$message <- sub("cliff_delta(): ",
                       paste0("randomization_test(): for the Cliff's delta ",
                              "rows (computed by cliff_delta()), "),
                       w$message, fixed = TRUE)
      warning(w)
      invokeRestart("muffleWarning")
    })
  cd_hat   <- cd$value[cd$term == "cliff_delta"]
  cd_lower <- cd$value[cd$term == "ci_lower"]
  cd_upper <- cd$value[cd$term == "ci_upper"]

  out <- data.frame(
    term = c("mean_difference", "statistic", "p_value", "p_value_se",
             "shift_ci_lower", "shift_ci_upper",
             "normal_theory_ci_lower", "normal_theory_ci_upper",
             "smd", "smd_ci_lower", "smd_ci_upper",
             "cles", "cles_ci_lower", "cles_ci_upper",
             "cliff_delta", "cliff_delta_ci_lower", "cliff_delta_ci_upper",
             "n_1", "n_2", "N", "n_evaluated", "exact"),
    value = c(mean_difference, T_obs, p_value, p_se,
              shift_lower, shift_upper,
              nt_lower, nt_upper,
              d, d_lower, d_upper,
              cl_hat, cl_lower, cl_upper,
              cd_hat, cd_lower, cd_upper,
              n_1, n_2, N, n_eval, as.integer(use_exact)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  out <- .as_dmar_tbl(out, conf_level = conf_level, p_terms = "p_value")
  attr(out, "statistic_name") <- if (statistic == "mean")
    "difference in means" else "studentized (Welch) t"
  attr(out, "method") <- if (use_exact) "exact enumeration" else "Monte Carlo"
  attr(out, "alternative")            <- alternative
  attr(out, "group_labels")           <- lev
  attr(out, "response_name")          <- response_name
  attr(out, "group_name")             <- group_name
  attr(out, "seed")                   <- seed
  attr(out, "observed_statistic")     <- T_obs
  attr(out, "reference_distribution") <- T_null
  attr(out, "reference_extreme")      <- flags
  out
}

# ---- internals ----
#
# The refusal of statistic = "t" when both observed groups are constant at
# different values. The message states why the statistic has no value,
# where the observed data lie in the reference set, and the exact p-value
# of the difference in means for the alternative of the call. With every
# score of one group equal to a and every score of the other equal to b,
# the observed assignment is the only one that reaches the observed
# difference in its own direction (1 / choose(N, n_1)); the assignment
# that exchanges the groups reaches it in absolute value, and exists only
# when n_1 = n_2 (2 / choose(N, n_1), otherwise 1 / choose(N, n_1)); and
# every assignment lies at or beyond it in the other direction (1).

.rt_refuse_constant_groups <- function(y_1, y_2, alternative, n_possible,
                                       use_exact, forced, n_resamples) {
  n_1 <- length(y_1)
  n_2 <- length(y_2)
  difference <- mean(y_1) - mean(y_2)
  toward <- if (difference > 0) "greater" else "less"
  numerator <- if (alternative == toward) 1 else
    if (alternative == "two_sided") (if (n_1 == n_2) 2 else 1) else n_possible
  p_value_mean <- numerator / n_possible
  reading <- if (numerator < n_possible)
    sprintf(paste0("the smallest p-value this alternative can give with ",
                   "groups of %d and %d"), n_1, n_2)
  else
    sprintf(paste0("because every reassignment has a mean difference at ",
                   "or %s the observed one"),
            if (difference > 0) "below" else "above")
  constant_in <- if (n_1 == n_2)
    paste0("only the observed one and the one that exchanges the two ",
           "groups leave both groups constant")
  else "only the observed one leaves both groups constant"
  monte_carlo <- if (!use_exact && numerator < n_possible)
    .rt_monte_carlo_text(sprintf("%d scores", n_1 + n_2), forced,
                         n_resamples)
  else ""
  .dmar_refuse(c(
    "randomization_test(): statistic = \"t\" is not defined for these ",
    "data. Every score in the first group is ", .rt_number_text(mean(y_1)),
    " and every score in the second is ", .rt_number_text(mean(y_2)),
    ", so neither group has any spread, the separate-variances standard ",
    "error of the mean difference is 0, and the studentized statistic has ",
    "a zero denominator. The observed data lie at the boundary of the ",
    "reference set: of the ", .rt_count_text(n_possible),
    " reassignments of the ", n_1 + n_2, " scores, ", constant_in,
    ", and every other reassignment has a positive standard error and a ",
    "finite statistic. The statistic therefore has no value for the data ",
    "that were observed, and no limiting value (Inf) is reported in its ",
    "place. This is a statement about the statistic and not about the ",
    "evidence. The difference in means is defined for every reassignment, ",
    "and statistic = \"mean\" gives the mean difference ",
    .rt_number_text(difference), " with the exact p-value ",
    .rt_p_text(numerator, n_possible), " (", alternative, "), ", reading,
    ".", monte_carlo, " Use statistic = \"mean\"."),
    "studentized_statistic_undefined",
    reason = "constant_groups", alternative = alternative,
    p_value_mean = p_value_mean, n_1 = n_1, n_2 = n_2)
}
