#' Factor Scores With Standard Errors, Intervals, and Determinacy
#'
#' Computes, for every row given to a fitted \pkg{lavaan} model and every
#' latent variable, both standard factor scores: the regression score
#' (Thomson, 1934; Thurstone, 1935) and the Bartlett score (Bartlett,
#' 1937, 1938). Each comes with a per-case standard error and confidence
#' interval, beside the number of the factor's items the case answered,
#' the case's missing data pattern, and the determinacy of each score for
#' that pattern (the correlation of the score with the factor it stands
#' for). For ordered items the two scores are their categorical
#' counterparts in \pkg{lavaan}, the empirical Bayes modal (EBM) score and
#' the maximum likelihood (ML) score.
#'
#' Both scores are returned because neither is better in general; which
#' one belongs in a later analysis depends on the role the factor plays
#' there (see \emph{Which Score Belongs Where}). For relations among
#' factors, the scores are rarely the right tool at all: a latent variable
#' model, estimated in one step or by structural after measurement
#' (\code{lavaan::sam()}; Rosseel & Loh, 2024), estimates them without the
#' bias a regression among scores carries. The scores are for what only
#' scores can do: describe and rank individual cases, carry a factor into
#' an analysis outside the model, or check the model's assumptions case by
#' case. Requires \pkg{lavaan}.
#'
#' @param fit A fitted \pkg{lavaan} model (from \code{lavaan::cfa()},
#'   \code{lavaan::sem()}, or \code{\link{fa_confirmatory}(output = "fit")}) estimated
#'   from raw data. A fit to summary statistics (it holds no rows) and a
#'   fit that did not converge are refused, as are a two-level fit, a fit
#'   of continuous items that keeps incomplete cases without FIML
#'   (two-stage maximum likelihood, pairwise deletion or the available
#'   cases) and is conditional
#'   on exogenous covariates, a fit with regressions among the variables
#'   \pkg{lavaan} holds as latent in which the disturbance of an observed
#'   variable regressed on others covaries with a factor (an observed
#'   mediator, for one), a
#'   fit conditional on exogenous covariates with regressions and an
#'   observed variable that \pkg{lavaan} holds among the latent variables
#'   but regresses on nothing, and, with ordered items, a fit with a direct
#'   effect of an exogenous covariate on an observed variable that also
#'   depends on a factor directly, one with a regression on an observed
#'   variable that is not an exogenous covariate when the variable predicted
#'   also depends on a factor directly, the predictor is ordered, or, in a
#'   row the fit
#'   keeps, the predictor is missing where the variable predicted is
#'   observed and depends on a factor through no variable the row has
#'   observed, one with an
#'   observed variable whose disturbance covaries with a factor and whose
#'   row \pkg{lavaan} keeps among the loadings, one with exogenous
#'   covariates that is not conditional on them, and any fit under a
#'   \pkg{lavaan} older than 0.7-3
#'   (see \emph{Details}). A fit with ordered items and residual
#'   covariances is scored with a warning.
#'   For continuous unconditional fits, a disturbance related only to an
#'   observed conditioning variable leaves the factor's regression
#'   conditional mean valid. Its scores and intervals are retained;
#'   Bartlett scores, SEs, limits, and determinacy are \code{NA} with the
#'   reason when their stronger conditional-unbiasedness requirement fails.
#' @param data Optional but recommended: the data frame given to the
#'   model, with the same rows in the same order (extra columns are
#'   ignored). It is checked against the data the model holds (every
#'   observed variable present, the same number of rows, and the same
#'   values in the rows the model used), and a frame that differs is
#'   refused. With it, the table has one row for every row of \code{data},
#'   including rows the fit did not use (deleted listwise, or with a
#'   missing group), and their reason and items answered are filled in;
#'   for a multiple group fit it must hold the group variable. Without
#'   it, the table runs from row 1 to the last row the fit holds, so rows
#'   at the end of the data that the fit did not use are not represented,
#'   and the reason given a row the fit did not use cannot say which of its
#'   values was missing.
#'   The scores are always those of the rows the fit was estimated from;
#'   scoring new observations is not what this function does.
#' @param conf_level Confidence level for the per-case intervals; defaults
#'   to 0.95.
#'
#' @details
#' \strong{The Two Scores.} For the observed items \eqn{o} of a case, with
#' \eqn{\Lambda_o} their rows of the loading matrix, \eqn{\Theta_{oo}}
#' their residual covariance matrix, \eqn{\Phi} the model implied
#' covariance matrix of the factors, \eqn{\Sigma_{oo} = \Lambda_o \Phi
#' \Lambda_o' + \Theta_{oo}}, and \eqn{\mu_o} the model implied
#' item means,
#' \deqn{\hat{\eta}_R = \mathrm{E}(\eta) + \Phi \Lambda_o'
#'   \Sigma_{oo}^{-1} (y_o - \mu_o)}
#' is the regression score, the mean of the factor's distribution given
#' the observed items under the normal factor model, and
#' \deqn{\hat{\eta}_B = \mathrm{E}(\eta) + (\Lambda_o'
#'   \Theta_{oo}^{-1} \Lambda_o)^{-1} \Lambda_o' \Theta_{oo}^{-1}
#'   (y_o - \mu_o)}
#' is the Bartlett score, the generalized least squares (and, under
#' normality, maximum likelihood) estimate of the case's factor value
#' treated as a fixed quantity. The regression score is shrunk toward the
#' factor mean, by more when the items carry less information, and so has
#' the smaller mean squared error; the Bartlett score is not shrunk, so
#' its expected value given the case's factor value is that value
#' (conditionally unbiased; under missing data this needs a condition on
#' why values are missing, see \emph{Missing Data}), at the price of a
#' larger variance. Their
#' standard errors differ in kind: that of the regression score is the
#' posterior standard deviation, the square root of the diagonal of
#' \eqn{\Phi - \Phi \Lambda_o' \Sigma_{oo}^{-1} \Lambda_o \Phi},
#' and that of the Bartlett score is its sampling standard error given the
#' case's factor value, the square root of the diagonal of
#' \eqn{(\Lambda_o' \Theta_{oo}^{-1} \Lambda_o)^{-1}}. The scores
#' and their standard errors are those of \code{lavaan::lavPredict()} with
#' \code{se = "standard"}, aligned to the rows of the data (with auxiliary
#' variables, the same formulas applied to the indicators alone; see
#' below).
#'
#' The Bartlett score of a factor exists only when the observed items
#' identify the factor's value, that is, when the factor's unit vector
#' lies in the row space of \eqn{\Lambda_o}, so that its value is a
#' function of the items' expected values. With simple structure that
#' holds for every factor with an observed item. With items shared among
#' factors (cross-loadings, a bifactor model) a pattern can fail it: a case
#' that answered one item loading on two factors, and nothing else of
#' either, has items that do not separate the two. The matrix
#' \eqn{\Lambda_o' \Theta_{oo}^{-1} \Lambda_o} is then singular, and
#' \pkg{lavaan} falls back on a generalized inverse, whose minimum norm
#' solution is not an estimate of either factor's value and whose standard
#' error is not the error of anything; such a factor's Bartlett score,
#' standard error, and determinacy are \code{NA} here with the reason. A
#' factor the pattern does identify keeps its score, which the generalized
#' inverse gives uniquely.
#'
#' \strong{Which Score Belongs Where.} Using scores in place of factors in
#' a later regression (factor score regression) with the same type of
#' score on both sides gives inconsistent estimates of the structural
#' slopes. Skrondal and Laake (2001) showed that the revised scheme, the
#' regression score for an explanatory factor and the Bartlett score for
#' a response factor, gives consistent estimates of every parameter of a
#' latent regression model, provided the explanatory and the response
#' factors are fitted as separate measurement models; scores from one
#' joint model of both are not covered by that result and should not feed
#' a regression among the same factors, and fitting each
#' factor alone (factorwise) is consistent only in special cases such as
#' the independent clusters model with one explanatory factor. Even then
#' the standard errors the second step reports do not account for the
#' scores having been estimated (corrected ones may require resampling),
#' no conventional test of fit is available, and factor score regression
#' works badly for nonlinear models such as polynomial regression.
#' Devlieger, Mayer, and Rosseel (2016) compared regression and Bartlett
#' factor score regression, the bias avoiding method of Skrondal and
#' Laake, and the bias correcting method of Croon (2002). The bias
#' avoiding method was unbiased for unstandardized but not for
#' standardized coefficients and had the largest standard error bias of
#' the four, and Croon's method, with the standard error they developed
#' for it, was the only suitable alternative to structural equation
#' modeling, comparable in bias, efficiency, mean squared error, power,
#' and Type I error rate, though with more standard error bias. The
#' differences among the methods disappeared when factor score
#' indeterminacy was low, and they advise assessing it first, which the
#' determinacy columns here allow. Croon's correction is local
#' structural after measurement with the maximum likelihood mapping, whose
#' matrix is the Bartlett score matrix (Rosseel & Loh, 2024), and Rosseel
#' and Loh recommend local structural after measurement over one-step
#' estimation when the sample is small to medium or the items have low to
#' moderate reliability, with as many separate measurement blocks as the
#' model permits. It is available as \code{lavaan::sam()}. Rosseel and
#' Loh developed the framework for complete continuous data, so with
#' missing data or ordered items its two-step standard errors deserve a
#' check against a bootstrap.
#'
#' \strong{Missing Data.} Under full information maximum likelihood
#' (\code{missing = "ml"}), \pkg{lavaan} scores each case from its own
#' observed items, pattern by pattern, which is what \eqn{o} denotes
#' above. For continuous items, the regression score so computed is the
#' posterior weighted maximum likelihood estimate that Estabrook and Neale
#' (2013) call full ML, and the Bartlett score so computed is their
#' unweighted ML estimate. In their simulation (one factor), the two were
#' equally reliable with complete data, and with missing values the
#' posterior weighted score was the more reliable (with five items and
#' values missing completely at random, a mean reliability of .780 against
#' .767, and .752 for the scaled person mean), which they attribute to its
#' empirical Bayes component. A case that answered none of a factor's
#' items has no Bartlett score for that factor; its regression score is
#' kept, but it comes only from the case's other observed variables,
#' through the model (the items of correlated factors, a covariate that
#' predicts the factor, an outcome the factor predicts), it is the
#' factor's mean when none of them bears on the factor, and its
#' determinacy says how little that is. A case with no observed variable
#' has no score of either kind. Under listwise deletion the rows with a
#' missing value were not used by the fit and are not scored; for
#' continuous items, a fit with \code{missing = "ml"} scores them, unless
#' the missing value is an exogenous covariate's. With \code{fixed.x =
#' TRUE}, \pkg{lavaan}'s default, every fit deletes a row whose covariate
#' is missing, FIML with \code{missing = "ml"} included, except one with
#' \code{missing = "ml.x"} (which a fit conditional on the covariates
#' does not allow), and the reason column names that remedy; with ordered
#' items no fit that \pkg{lavaan} can score keeps such a row.
#' \pkg{lavaan} has no FIML for ordered items: fitted by WLSMV, they take
#' their parameters from pairwise or listwise deletion, which assumes the
#' values are missing completely at random, and under pairwise deletion
#' (\code{missing = "pairwise"}) the EBM and ML scores use each case's
#' observed items, so that fit scores the rows listwise deletion leaves
#' out.
#'
#' A fit of continuous items that keeps the incomplete cases without FIML
#' is scored the same way, each case from its own observed items, with the
#' fit's estimates in the formulas above: a fit by two-stage maximum
#' likelihood (\pkg{lavaan}'s \code{missing = "two.stage"} or
#' \code{"robust.two.stage"}, as \code{\link{fa_confirmatory}} fits with
#' \code{missing = "two_stage"}) or by pairwise deletion (\code{missing =
#' "pairwise"}, or \code{"available.cases"}, which for continuous items
#' gives \pkg{lavaan} the same pairwise complete moments and estimates).
#' The two-stage estimates are consistent when the values
#' are missing at random, as the FIML estimates are (Savalei & Bentler,
#' 2009); the estimates from pairwise deletion require the values to be
#' missing completely at random. The scores treat the estimates as known,
#' as they treat the FIML estimates, so the scores, standard errors, and
#' determinacies of such a fit mean what they mean under FIML and differ
#' from FIML's only as the estimates differ. A fit without a mean
#' structure (the default for one group under pairwise deletion) centers
#' its observed variables at their sample means and its factors at the
#' means \pkg{lavaan} gives them (\code{lavInspect(fit, "mean.lv")}): zero,
#' or for a factor a covariate predicts, the prediction at the
#' covariate's mean, as a FIML fit has it. \code{lavaan::lavPredict()}
#' scores such a fit as complete data, leaving every case with a missing
#' value without a score, so these scores are computed here, by the closed
#' forms \code{lavPredict()} applies pattern by pattern under FIML.
#' Auxiliary variables of the first stage of a two-stage fit
#' (\pkg{lavaan}'s \code{aux}) are not in the model and do not enter the
#' scores. Such a fit conditional on exogenous covariates
#' (\code{conditional.x = TRUE}) is refused.
#'
#' Missing values affect the two scores differently. Conditional
#' unbiasedness of the Bartlett score and the coverage of its Wald
#' interval are properties of its sampling distribution given the case's
#' factor value, and with missing values they hold when, given the
#' factors, whether an item is observed does not depend on the residuals
#' of the items scored: values missing completely at random, or missing
#' because of the factors themselves or of variables related to the items
#' only through the factors. Rubin (1976) showed that sampling
#' distribution inferences may ignore the process that causes missing
#' data when the missing values are missing at random and the observed
#' values are observed at random, whereas likelihood and Bayesian
#' inferences need only the first (with the parameters of the process
#' distinct from those of the data). Missingness that depends on the
#' observed values of the factor's own items, the kind of missing at
#' random that FIML is commonly used for, fails the second: it selects
#' each pattern on those items' residuals, so the Bartlett score of a
#' pattern is biased given the factor value, and its interval covers the
#' factor value less often than stated in that pattern. In a one-factor
#' model of four items with loadings of .7 and 5,000 cases, the last three
#' items missing whenever the first fell below -0.3, the Bartlett scores
#' of the cases left with one item averaged 0.82 below their factor
#' values and their 95\% intervals covered 92.7\% of them; with as many
#' values missing completely at random, -0.01 and 95.3\%. The regression
#' score is a posterior quantity, so missingness that depends only on
#' what the score conditions on leaves its posterior interval with its
#' average coverage (94.6\% in the first case). Missingness that runs
#' through an auxiliary variable the scores set aside (see
#' \emph{Auxiliary Variables and Covariates}) is not of that kind for the
#' regression score, since the pattern then carries information about the
#' factor beyond the observed items, and the regression score of a pattern
#' is biased. The Bartlett score is not, but only when the auxiliary is
#' related to the items through the factors alone (independent of them
#' given the factors); one that shares variance with an item beyond the
#' factor, as an earlier wave of the same item does, selects the pattern
#' on that item's residual, which fails the condition above. In the same
#' design with the last three items missing whenever an auxiliary fell
#' below -0.3, an auxiliary that correlated .4 with the factor and with
#' nothing else left the cases with one item Bartlett scores 0.05 below
#' their factor values on average (within two Monte Carlo standard errors
#' of zero), with intervals covering 95.6\%, and regression scores 0.17
#' above them; one that correlated .4 with the factor and .5 with the first
#' item's residual left them Bartlett scores 0.56 below their factor
#' values, with intervals covering 93.1\%. A score whose
#' error has mean zero given the factor leaves the slopes of a regression
#' in which it is the response unbiased, and that is the Bartlett score's
#' place in the revised scheme of factor score regression above, so
#' missingness that depends on the factor's own observed items weakens it
#' there. The reliabilities of Estabrook and Neale (2013) quoted above are
#' for values missing completely at random.
#'
#' \strong{Determinacy.} The determinacy of a score is its correlation
#' with the factor, \eqn{[A \Lambda_o \Phi]_{kk} / \sqrt{[A \Sigma_{oo}
#' A']_{kk} \Phi_{kk}}} for the weight matrix \eqn{A} of the
#' score, the formula \pkg{lavaan} uses for \code{fs.determinacy} and
#' \code{fs.determinacy.bartlett}, evaluated here for each missing data
#' pattern. Its square is the reliability of the score as a measure of the
#' factor. The regression score has the larger determinacy, since it is
#' the linear function of the observed items that correlates most with
#' the factor; for a one-factor model the two are equal, and with
#' uncorrelated residuals their square is coefficient \emph{H} (Hancock &
#' Mueller, 2001), the reliability of the optimally weighted composite of
#' the items. Per-pattern determinacy
#' shows the cost of each missing item. It is computed for continuous
#' items; with ordered items each case's standard error is given instead.
#' Beside the regression score's determinacy \eqn{\rho}, the display and
#' the \code{determinacy} and \code{patterns} attributes report
#' \eqn{2\rho^2 - 1} (\code{min_score_correlation}), the lowest
#' correlation that two sets of scores, each as consistent with the model
#' as the other, can have (Guttman, 1955, as Grice, 2001, p. 435, presents
#' it). It is zero when \eqn{\rho^2 = .5}, where two such sets can be
#' uncorrelated, and negative below that, so it says how far the scores of
#' different cases can be trusted to order the cases on the factor. The
#' package reports the quantity and draws no line.
#'
#' \strong{Intervals.} The interval around a regression score is a
#' posterior (empirical Bayes) interval: its coverage holds on average
#' over cases, not for each case's factor value. The interval around a
#' Bartlett score is a Wald interval for the case's factor value. Both
#' treat the model's estimated parameters as known, so they do not carry
#' the parameters' sampling variability. For ordered items the EBM
#' interval approximates the posterior around its mode, and the ML
#' interval is a Wald interval. The ML score of a factor exists only when
#' the likelihood of the case's responses has a finite maximum in it. A
#' case whose observed items of a factor all sit at the same end of their
#' scales (all in the lowest category, or all in the highest, for positive
#' loadings) has none, since the likelihood keeps rising as the factor's
#' value moves out. With items shared among factors (cross-loadings, a
#' bifactor model) the likelihood can also rise without bound along a
#' direction that moves several factors together, or stay flat along one
#' that the observed items cannot separate. The check is exact: a factor
#' has a finite and unique ML score when no direction along which no
#' observed item's probability falls moves it, which a nonnegative least
#' squares problem decides (by Farkas' lemma), and with simple structure it
#' reduces to the rule for the items at an end. It reads each item's
#' loadings as \pkg{lavaan}'s likelihood does, including those of an
#' observed variable regressed on the factors (an outcome the factor
#' predicts, or an item with a direct effect of another factor), which
#' \pkg{lavaan} holds among the regressions. A factor that fails it has
#' its ML score and interval \code{NA} with the reason. A factor that
#' passes it keeps its ML score, also when it shares observed items with a
#' factor whose likelihood rises without bound: as that factor moves out,
#' the shared items' probabilities tend to 1 whatever this factor's value,
#' so its maximum is that of the case's other items, which \pkg{lavaan}'s
#' optimizer reaches within its tolerance (within 1.2e-4 of a direct
#' maximization of the likelihood for the factor of a cross-loading item
#' in the tests, and within 1.1\% of the standard error in a bifactor
#' model of 12 items).
#' When \pkg{lavaan}'s default optimizer stops on its iteration limit and
#' returns no score for the case, the score comes from \pkg{lavaan}'s
#' \code{nlminb} optimizer.
#' A higher order factor, which has no items of its own, has no ML score at
#' all: the likelihood of the responses depends on the factors it
#' measures, not on it, and \pkg{lavaan}'s value is the center of the small
#' ridge its optimizer adds. Its EBM score, which the factor's place in the
#' model informs, is kept.
#'
#' \strong{Ordered Items Need lavaan 0.7-3.} The probability of an ordered
#' response compares each threshold with the mean of the item's latent
#' response, its intercept plus its loadings times the factor values. The
#' EBM and ML scores of \pkg{lavaan} 0.7-2 counted the intercept twice in
#' that comparison, so they were wrong in any model whose ordered
#' indicators have intercepts other than zero (a multiple group model with
#' equal thresholds, for one); \pkg{lavaan} 0.7-3 corrected it, and also
#' changed how the scores treat the scale factors of the delta
#' parameterization. A fit with ordered items is therefore refused under an
#' older \pkg{lavaan}, with the remedy, rather than scored with numbers
#' that may be wrong.
#'
#' \strong{Multiple Groups.} The scores use each group's own parameters,
#' and the \code{group} column names the group. Scores from different
#' groups are on a common metric only to the extent the measurement model
#' is invariant across the groups (equal loadings for comparing relations
#' that involve the scores, and equal intercepts as well for comparing
#' their levels); see
#' \code{\link{measurement_invariance}}. Rows whose group is missing are
#' not used by the fit; with \code{data} they are kept, without a score.
#'
#' \strong{Auxiliary Variables and Covariates.} An observed variable that
#' indicates no factor and enters the model only through covariances with
#' the indicators is an auxiliary variable carried as a saturated
#' correlate (Graham, 2003), as the reliability functions of this package
#' carry them. It helps the estimates under missing data, but under that
#' parameterization it is uncorrelated with every factor and correlated
#' with the item residuals, so whatever it shares with the factor is
#' carried as residual covariance, and a score that conditions on it
#' subtracts true score variance. When such variables are present the
#' scores are computed from the indicators alone, with the fit's
#' estimates, the \code{scored_from} attribute is \code{"indicators"}, and
#' a message says so. A higher order factor, which has no items of its
#' own, then has no Bartlett score, nor, when the data have missing values,
#' in the scores computed here for a two-stage fit or one by pairwise
#' deletion (\code{missing = "pairwise"} or \code{"available.cases"}), nor
#' in a fit in which \pkg{lavaan} gives an observed variable a row of zeros
#' among the loadings (see
#' \emph{Observed Variables in Regressions}): \pkg{lavaan} forms that
#' score through the loadings of the factors it measures, a construction
#' not repeated here, and the reason column says so; its regression score
#' is computed as for any other factor. A covariate that predicts a
#' factor (\code{f ~ x}) is different: the scores condition on it as the
#' model specifies (as on any observed variable that enters a regression
#' of the model), and Curran, Cole, Bauer, Hussong, and Gottfredson (2016)
#' found that including such covariates in the scoring model improved
#' score quality in nearly every condition of their simulation.
#'
#' \strong{Observed Variables in Regressions.} An observed variable regressed
#' on others, such as an item with a direct effect of a covariate
#' (\code{y1 ~ x}, the multiple indicators multiple causes model of
#' differential item functioning) or an outcome the factor predicts
#' (\code{z ~ f}), is held by \pkg{lavaan} as a latent variable measured by
#' itself without error, with its loadings among the regressions and its
#' residual variance among those of the latent variables. The scores, their
#' standard errors and determinacies, the check of the ML score, and the
#' check of the solution take its loadings and residual variance from there,
#' as \code{lavaan::lavPredict()} does, so an item with a direct effect
#' enters its factor's Bartlett score, which is then conditionally unbiased
#' given the factor and the covariate. An observed variable that \pkg{lavaan}
#' holds the same way but regresses on nothing, such as a criterion that
#' covaries with a factor (\code{f ~~ z}) or a covariate of a factor that
#' also covaries with another variable (so that \pkg{lavaan} does not count
#' it among the exogenous covariates), is conditioned on, as an exogenous
#' covariate is: the regression score is the factor's mean given it and the
#' items, and the Bartlett score is conditionally unbiased given the factors
#' and it. In a fit conditional on the covariates (\code{conditional.x =
#' TRUE}) without regressions among the variables \pkg{lavaan} holds as
#' latent (no \eqn{B} matrix), so is an observed variable regressed only on
#' exogenous covariates (\code{z ~ x}): \pkg{lavaan} then keeps the
#' variable's own row among the loadings, \code{lavPredict()}'s regression
#' score is the factor's mean given the items, the variable, and the
#' covariates, and its Bartlett score is conditionally unbiased given the
#' factors, the variable, and the covariates, also when the variable
#' covaries with a factor (\code{f ~ x}, \code{z ~ x}, \code{f ~~ z};
#' \code{lavaan::sem()} frees that covariance by default). In a model with
#' regressions among them, \code{lavPredict()} gives a variable regressed
#' on nothing a row of zeros among the loadings, which leaves its covariances
#' with the factors out of the scores (with \code{f ~ x} and \code{f ~~ z}
#' for three Holzinger and Swineford tests, its regression score missed the
#' conditional mean by up to 1.9), so the scores of such a fit are computed
#' here, with the formulas \code{lavPredict()} uses; one conditional on the
#' covariates (\code{conditional.x = TRUE}) is refused, with the remedy. In a
#' fit without a mean structure and without exogenous covariates,
#' \code{lavPredict()} counts the mean of a variable held this way twice, as
#' its intercept and as the mean of the latent variable that holds it, which
#' shifts by a constant every score whose weights involve the variable: the
#' regression score of every factor related to it and the Bartlett score of a
#' factor that predicts it (for a factor of three Holzinger and Swineford
#' spatial tests that predicts a fourth test and correlates .50 with a factor
#' of three speed tests, by -0.39 and -0.61 for the first factor's regression
#' and Bartlett scores and by -0.07 for the second factor's regression
#' score); the scores of such a fit are computed here, with the formulas
#' \code{lavPredict()} uses and the means the fit implies. A higher order
#' factor's Bartlett score, which \pkg{lavaan} forms through the loadings of
#' the factors it measures, is \pkg{lavaan}'s less the constant the double
#' count adds (-1.09 for a second order factor of three factors of three
#' Holzinger and Swineford tests each, when the spatial factor predicts the
#' paper form board test), so it equals that of the fit with a mean
#' structure. In a model with regressions among the variables \pkg{lavaan}
#' holds as latent, a variable regressed on others is represented by its
#' measurement equation, with its disturbance as the error, which is right
#' only when that disturbance is uncorrelated with the factors and with the
#' variables the scores condition on, so a fit in which such a variable
#' predicts a factor, directly or through others (an observed mediator), or
#' has a residual covariance with a factor is refused: \code{lavPredict()} then gives a Bartlett score
#' that is not conditionally unbiased and, when a factor is involved, a
#' regression score that is not the factor's conditional mean (with
#' \code{f ~ a} and \code{a ~ x} for three Holzinger and Swineford tests,
#' it missed the conditional mean by up to 1.08, with a correlation of
#' .92). Carried as a
#' single indicator latent variable (\code{A =~ a} with \code{a ~~ 0*a}),
#' such a variable is scored correctly: written so, \pkg{lavaan} fixes one
#' scale constraint, the loading or, with \code{std.lv = TRUE}, the latent
#' variable's variance, so the refitted model is the same model, with the
#' refused fit's log likelihood and degrees of freedom (Bollen, 1989,
#' Chapter 7). Written \code{A =~ 1*a} in a fit with \code{std.lv = TRUE},
#' it fixes both, which is a different model whose scores are not those of
#' the refused one. In a multiple group fit, \code{group.equal} holds
#' parameters equal by kind, and the variable's intercept, residual
#' variance, and covariances become the latent variable's, so the refit is
#' the same model only when it carries the refused fit's constraints
#' across the groups: with equal intercepts, \pkg{lavaan} frees the latent
#' variable's mean after the first group unless it is fixed at zero in
#' every group (\code{A ~ 0*1}; for two schools with equal loadings and
#' intercepts the refit without it had a log likelihood of -4246.42 on 15
#' degrees of freedom, against -4246.51 on 16), and the equalities
#' \code{group.equal} places on residual variances and covariances, or on
#' those of latent variables, must be given to the variable's parameters as
#' the refused fit had them and to no others. The message says which
#' apply.
#'
#' With ordered items, \pkg{lavaan}'s EBM and ML scores leave a direct
#' effect of an exogenous covariate on an observed variable out of each
#' case's likelihood, and they take an observed variable that is not an
#' exogenous covariate as zero where it predicts another. That leaves the
#' scores unchanged when the variable predicted depends on no factor
#' directly: its density given its predictors then does not involve the
#' factors, so the term left out changes each case's likelihood by a factor
#' constant in them. Such fits are scored: a covariate's effect on an
#' observed outcome unrelated to the factors (\code{w ~ x}), and a fully
#' mediated chain of continuous observed variables (\code{m ~ f} with
#' \code{w ~ m}, whose EBM and ML scores agree with a direct maximization
#' of the model's likelihood for four Holzinger and Swineford tests cut at
#' their quartiles). A fit in which the term left out involves the factors
#' is refused: a direct effect on a variable that also depends on a factor
#' directly (\code{y1 ~ x} for an item, whose thresholds are estimated
#' given the
#' covariate, so the scores are those of a case with the covariate at zero;
#' ML scores off by 0.11 to 0.44 for an item regressed on age), a
#' regression of such a variable on an observed one (\code{y4 ~ y3} for
#' local dependence between two items, ML scores off by up to 0.18 for the
#' four tests; or \code{w ~ m + f}), a regression on an ordered item, whose
#' latent response is not observed (\code{w ~ y3}), and a regression on an
#' observed variable that a row the fit keeps lacks while it has the one
#' predicted, where the missing variable depends on a factor through no
#' variable the row has observed (pairwise deletion, where the model's
#' likelihood of the variable predicted then runs through the missing one
#' to the factors). A variable the row has observed on the way cuts that
#' dependence, and the row is scored: with \code{m ~ f}, \code{w ~ m},
#' and \code{s ~ w}, a row lacking \code{w} but having \code{m} and
#' \code{s} is scored, while a row lacking \code{m} but having \code{w}
#' refuses the fit. A fit is refused, too, for an observed variable
#' \pkg{lavaan} holds among the
#' latent variables whose disturbance covaries with a factor, regressed on
#' nothing or, in a fit without regressions among the latent variables,
#' only on exogenous covariates (\code{f ~ x} with \code{w ~ x}, where
#' \code{lavaan::sem()} frees the covariance of \code{f} and \code{w} by
#' default): the EBM score takes the factor's prior from its model implied
#' mean given the covariates, which leaves the variable's value out, so it
#' is the posterior mode given the items and covariates alone; written as a
#' predictor of the factor (\code{f ~ x}), a variable regressed on nothing
#' is an exogenous covariate, which the EBM score conditions on. A fit
#' with ordered items and exogenous covariates that is not conditional on
#' them (\code{conditional.x = FALSE}), which \code{lavPredict()} cannot
#' score, is refused as well. The EBM and ML scores of \pkg{lavaan}
#' also take each case's
#' likelihood as the product of its observed variables' probabilities given
#' the factors, so they leave residual covariances among those variables out;
#' a fit with one is scored with a warning, since its scores are those of the
#' model without it (with a residual correlation of -0.12 between two of the
#' four tests, the ML scores of the 160 response patterns differed from the
#' scores that maximize the model's own likelihood by up to 0.11, 15\% of
#' their standard error, and by 0.03 at the median).
#'
#' An improper solution is scored with a warning: one with a negative
#' residual variance (a Heywood case, also of an observed variable
#' regressed on others), a negative factor variance, a covariance matrix
#' of the factors that is not positive definite (a factor correlation
#' beyond 1, for one), a residual covariance matrix of the continuous
#' variables scored from that is not positive definite, or, when none of
#' these holds, a residual covariance matrix of the structural part of the
#' model that is not, the conditions \pkg{lavaan}'s post-estimation check
#' flags. The regression weights use the model implied covariance matrix
#' and the Bartlett weights the loadings and the inverse of the residual
#' covariance matrix, so both scores inherit the solution. For ordered
#' items \pkg{lavaan} computes no EBM or ML score for a solution with a
#' negative residual variance, and the warning and the reason column say
#' so.
#'
#' @return A \code{data.frame} (class \code{dmar_factor_scores},
#'   \code{dmar_tbl}) with one row for each row of the data and each
#'   factor, ordered by factor and then row, so that
#'   \code{x$.derived_bartlett[x$factor == "f"]} lines up with the rows of
#'   \code{data}. Its columns are always the same: \code{row} (the row of
#'   the data), \code{group} (the group label; \code{NA} for a fit without
#'   groups), \code{factor}, \code{.derived_regression},
#'   \code{.derived_regression_se}, \code{.derived_regression_ci_lower},
#'   \code{.derived_regression_ci_upper}, \code{.derived_bartlett},
#'   \code{.derived_bartlett_se}, \code{.derived_bartlett_ci_lower},
#'   \code{.derived_bartlett_ci_upper}, \code{items_observed} and \code{items} (the
#'   factor's items the case answered, and its number of items; a higher
#'   order factor's items are those of the factors it measures, and a
#'   loading the model holds at zero, written \code{0*x} or \code{a*x}
#'   with \code{a == 0} directly or through a chain of equalities, makes
#'   no item: the scores and reasons are those of the model without it),
#'   \code{pattern} (the missing data pattern of the variables scored
#'   from, within the group, numbered from the most complete),
#'   \code{determinacy_regression}, \code{determinacy_bartlett} (for the
#'   case's pattern), and \code{reason}
#'   (empty, or why a value is missing or borrowed; reasons are separated
#'   by semicolons). The \code{.derived_regression} and
#'   \code{.derived_bartlett} columns hold the EBM and ML scores for
#'   ordered items.
#'
#'   The eight columns that hold the scores, their standard errors, and
#'   their confidence limits carry the prefix \code{.derived_}, so that a
#'   score joined to the observed data cannot be taken for an observed
#'   variable or collide with one (a data set may well have a column
#'   named \code{regression}). The \code{derived_columns} attribute
#'   names them by what each holds (\code{regression},
#'   \code{regression_se}, \code{regression_ci_lower},
#'   \code{regression_ci_upper}, and the same four for \code{bartlett}),
#'   so \code{x[[attr(x, "derived_columns")[["bartlett"]]]]} is the column
#'   of Bartlett scores. All are present under every supported model.
#'
#'   Attributes: \code{methods}
#'   (the \pkg{lavaan} methods behind the two columns), \code{conf_level},
#'   \code{estimator}, \code{missing}, \code{patterns} (one row per group,
#'   pattern, and factor: the number of cases, the items observed, the
#'   variables missing, both determinacies, \code{min_score_correlation},
#'   and both standard errors), \code{determinacy} (per group and factor,
#'   for complete data: both determinacies and \code{min_score_correlation}),
#'   \code{scored_from} (\code{"fit"} or \code{"indicators"}),
#'   \code{auxiliary} (the auxiliary variables set aside), \code{ngroups},
#'   and \code{lavaan_version}. \code{tidy()} returns the table as a plain
#'   data frame; \code{glance()} returns one row with the numbers of rows,
#'   rows scored, factors, groups, and patterns, the two methods, the
#'   estimator, the missing data treatment, \code{scored_from}, and
#'   \code{conf_level}.
#'
#' @references
#' Bartlett, M. S. (1937). The statistical conception of mental factors.
#'   \emph{British Journal of Psychology, 28}(1), 97--104.
#'
#' Bartlett, M. S. (1938). Methods of estimating mental factors.
#'   \emph{Nature, 141}(3570), 609--610.
#'
#' Bollen, K. A. (1989). \emph{Structural equations with latent
#'   variables}. Wiley.
#'
#' Croon, M. (2002). Using predicted latent scores in general latent
#'   structure models. In G. Marcoulides & I. Moustaki (Eds.),
#'   \emph{Latent variable and latent structure models} (pp. 195--223).
#'   Erlbaum.
#'
#' Curran, P. J., Cole, V., Bauer, D. J., Hussong, A. M., & Gottfredson,
#'   N. (2016). Improving factor score estimation through the use of
#'   observed background characteristics. \emph{Structural Equation
#'   Modeling, 23}(6), 827--844. \doi{10.1080/10705511.2016.1220839}
#'
#' Devlieger, I., Mayer, A., & Rosseel, Y. (2016). Hypothesis testing
#'   using factor score regression: A comparison of four methods.
#'   \emph{Educational and Psychological Measurement, 76}(5), 741--770.
#'   \doi{10.1177/0013164415607618}
#'
#' Estabrook, R., & Neale, M. (2013). A comparison of factor score
#'   estimation methods in the presence of missing data: Reliability and
#'   an application to nicotine dependence. \emph{Multivariate Behavioral
#'   Research, 48}(1), 1--27. \doi{10.1080/00273171.2012.730072}
#'
#' Graham, J. W. (2003). Adding missing-data-relevant variables to
#'   FIML-based structural equation models. \emph{Structural Equation
#'   Modeling, 10}(1), 80--100. \doi{10.1207/S15328007SEM1001_4}
#'
#' Grice, J. W. (2001). Computing and evaluating factor scores.
#'   \emph{Psychological Methods, 6}(4), 430--450.
#'   \doi{10.1037/1082-989X.6.4.430}
#'
#' Guttman, L. (1955). The determinacy of factor score matrices with
#'   implications for five other basic problems of common-factor theory.
#'   \emph{British Journal of Statistical Psychology, 8}(2), 65--81.
#'   \doi{10.1111/j.2044-8317.1955.tb00321.x}
#'
#' Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
#'   reliability within latent variable systems. In R. Cudeck, S. du
#'   Toit, & D. Sörbom (Eds.), \emph{Structural equation modeling:
#'   Present and future} (pp. 195--216). Scientific Software
#'   International.
#'
#' Rosseel, Y. (2012). lavaan: An R package for structural equation
#'   modeling. \emph{Journal of Statistical Software, 48}(2), 1--36.
#'   \doi{10.18637/jss.v048.i02}
#'
#' Rosseel, Y., & Loh, W. W. (2024). A structural after measurement
#'   approach to structural equation modeling. \emph{Psychological
#'   Methods, 29}(3), 561--588. \doi{10.1037/met0000503}
#'
#' Rubin, D. B. (1976). Inference and missing data. \emph{Biometrika,
#'   63}(3), 581--592. \doi{10.1093/biomet/63.3.581}
#'
#' Savalei, V., & Bentler, P. M. (2009). A two-stage approach to missing
#'   data: Theory and application to auxiliary variables. \emph{Structural
#'   Equation Modeling, 16}, 477--497. \doi{10.1080/10705510903008238}
#'
#' Skrondal, A., & Laake, P. (2001). Regression among factor scores.
#'   \emph{Psychometrika, 66}(4), 563--575. \doi{10.1007/BF02296196}
#'
#' Thomson, G. H. (1934). The meaning of \emph{i} in the estimate of
#'   \emph{g}. \emph{British Journal of Psychology, 25}, 92--99.
#'
#' Thurstone, L. L. (1935). \emph{The vectors of mind}. University of
#'   Chicago Press.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{fa_confirmatory}} to fit the measurement model;
#'   \code{\link{reliability_H}} for coefficient \emph{H};
#'   \code{\link{reliability_omega}} for the reliability of the
#'   unit-weighted composite; \code{\link{measurement_invariance}} before
#'   comparing scores across groups; \code{lavaan::lavPredict()},
#'   \code{lavaan::sam()}.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' data(holzinger_swineford)
#'
#' # Two later spatial tests (t25, t26) were given to the Grant-White
#' # pupils only, so under full information maximum likelihood each
#' # Pasteur pupil (rows 1 to 156) is scored from four spatial tests and
#' # each Grant-White pupil from six.
#' model <- "
#'   spatial =~ t1_visual_perception + t2_cubes + t3_paper_form_board +
#'              t4_lozenges + t25_paper_form_board_r + t26_flags
#'   verbal  =~ t5_general_information + t6_paragraph_comprehension +
#'              t7_sentence + t8_word_classification + t9_word_meaning"
#' fit <- lavaan::cfa(model, data = holzinger_swineford, std.lv = TRUE,
#'                    missing = "ml")
#' fs <- factor_scores(fit, data = holzinger_swineford)
#' fs
#'
#' # One pupil from each school: the spatial scores of the Pasteur pupil
#' # rest on four tests and have the wider intervals.
#' fs[fs$row %in% c(1, 157),
#'    c("row", "factor", ".derived_regression",
#'      ".derived_regression_ci_lower", ".derived_regression_ci_upper",
#'      ".derived_bartlett", ".derived_bartlett_ci_lower",
#'      ".derived_bartlett_ci_upper", "items_observed", "items")]
#'
#' # Determinacy and standard errors for each missing data pattern.
#' attr(fs, "patterns")
#'
#' # The table is ordered by factor and then row, so one factor's column
#' # lines up with the data.
#' spatial_bartlett <- fs$.derived_bartlett[fs$factor == "spatial"]
#' length(spatial_bartlett) == nrow(holzinger_swineford)
#'
#' # When scores must stand in for factors in a regression, Skrondal and
#' # Laake (2001) pair the regression score of the predictor with the
#' # Bartlett score of the outcome, each from its own measurement model.
#' # With both factor variances fixed at 1 the slope estimates the latent
#' # correlation; pairing the regression scores attenuates it. The lm()
#' # standard errors treat the scores as error free and are not valid.
#' fit_v <- lavaan::cfa("verbal =~ t5_general_information +
#'   t6_paragraph_comprehension + t7_sentence + t8_word_classification +
#'   t9_word_meaning", data = holzinger_swineford, std.lv = TRUE,
#'   missing = "ml")
#' fit_s <- lavaan::cfa("spatial =~ t1_visual_perception + t2_cubes +
#'   t3_paper_form_board + t4_lozenges + t25_paper_form_board_r +
#'   t26_flags", data = holzinger_swineford, std.lv = TRUE,
#'   missing = "ml")
#' fs_v <- factor_scores(fit_v)
#' fs_s <- factor_scores(fit_s)
#' c(revised = unname(coef(lm(fs_s$.derived_bartlett ~
#'                               fs_v$.derived_regression))[2]),
#'   regression_both = unname(coef(lm(fs_s$.derived_regression ~
#'                                       fs_v$.derived_regression))[2]),
#'   latent = lavaan::lavInspect(fit, "cor.lv")["spatial", "verbal"])
#'
#' # The broom verbs: the table as a plain data frame, and one summary row.
#' head(generics::tidy(fs))
#' generics::glance(fs)
#'
#' @export
factor_scores <- function(fit, data = NULL, conf_level = 0.95) {
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  }
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    .dmar_refuse(c(
      "Package 'lavaan' is required. Install it with ",
      "install.packages(\"lavaan\")."), "package_required", package = "lavaan")
  }
  .fs_check_fit(fit)

  opts <- lavaan::lavInspect(fit, "options")
  lv <- lavaan::lavNames(fit, "lv")
  ov <- lavaan::lavNames(fit, "ov")
  ordered <- lavaan::lavNames(fit, "ov.ord")
  categorical <- length(ordered) > 0L
  # lavaan 0.7-2's EBM and ML scores counted the intercepts of ordered
  # indicators twice (see .fs_require_lavaan_ordered()).
  if (categorical) .fs_require_lavaan_ordered()
  conditional_x <- isTRUE(opts$conditional.x)
  ngroups <- lavaan::lavInspect(fit, "ngroups")
  group_var <- as.character(lavaan::lavInspect(fit, "group"))
  glab <- if (length(group_var) == 1L && nzchar(group_var)) {
    as.character(lavaan::lavInspect(fit, "group.label"))
  } else {
    rep(NA_character_, ngroups)
  }
  index <- lapply(.fs_as_list(lavaan::lavInspect(fit, "case.idx")),
                  as.integer)
  model_data <- lapply(.fs_as_list(lavaan::lavInspect(fit, "data")),
                       as.matrix)
  # Each group's model matrices, with every entry the model holds at zero
  # set to exactly zero (.fs_zero_held()), and with lavaan's handling of the
  # dummy latent variables that carry observed variables in regressions
  # (.fs_dummies()); lavPredict() is given a copy of the fit with the same
  # exact zeros (`fit_lp`).
  ov_x <- lavaan::lavNames(fit, "ov.x")
  held <- .fs_held_cells(fit)
  est <- .fs_zero_held(.fs_as_list(lavaan::lavInspect(fit, "est"),
                                   inner = "lambda"), held)
  est <- lapply(est, .fs_dummies, lv = lv, ov_x = ov_x)
  fit_lp <- .fs_zero_held_fit(fit, held)
  methods <- if (categorical) {
    c(regression = "EBM", bartlett = "ML")
  } else {
    c(regression = "regression", bartlett = "Bartlett")
  }

  aux <- .fs_auxiliary(fit)
  if (length(aux) && categorical) {
    .dmar_refuse_argument(c(
      "Auxiliary variables carried as saturated correlates (",
      paste(aux, collapse = ", "), ") are a full information maximum ",
      "likelihood device for continuous items; with ordered items, ",
      "refit without them to score."), arg = NULL)
  }
  if (length(aux) && conditional_x) {
    .dmar_refuse_argument(c(
      "Auxiliary variables carried as saturated correlates (",
      paste(aux, collapse = ", "), ") cannot be set aside in a fit ",
      "conditional on exogenous covariates (conditional.x = TRUE); ",
      "refit with conditional.x = FALSE to score."), arg = NULL)
  }
  # Fits whose observed variables in regressions lavaan scores wrongly
  # (.fs_check_dummies()).
  bartlett_issues <- .fs_check_dummies(est, lv, glab, categorical, conditional_x, ov_x, aux,
                    ordered = ordered, data = model_data,
                    group_equal = as.character(opts$group.equal))
  # lavPredict() scores each case from its own observed items, pattern by
  # pattern, only under FIML (missing = "ml" or "ml.x"). A fit of
  # continuous items that keeps incomplete cases otherwise, by two-stage
  # maximum likelihood (missing = "two.stage" or "robust.two.stage"),
  # pairwise deletion (missing = "pairwise"), or the available cases, it
  # scores as complete data, so each case with a missing value gets no
  # score. Those scores are computed here instead, case by case from the
  # observed items with the fit's estimates, by the closed forms
  # lavPredict() applies pattern by pattern under FIML
  # (.fs_pattern_info()), as for a fit with auxiliary variables. Ordered
  # items are scored case by case by lavaan under any treatment of missing
  # values.
  per_case <- !categorical && !isTRUE(opts$missing %in% c("ml", "ml.x")) &&
    any(vapply(model_data, anyNA, logical(1)))
  fit_desc <- if (isTRUE(opts$missing %in% c("two.stage",
                                             "robust.two.stage"))) {
    c("A two-stage fit", "two-stage estimates")
  } else if (identical(opts$missing, "pairwise")) {
    c("A fit by pairwise deletion", "estimates from pairwise deletion")
  } else if (identical(opts$missing, "available.cases")) {
    c("A fit to the available cases", "estimates from the available cases")
  } else {
    c(sprintf("A fit with missing = \"%s\"", opts$missing),
      "fit's estimates")
  }
  if (per_case && conditional_x) {
    .dmar_refuse_argument(c(
      fit_desc[1L], " conditional on exogenous covariates ",
      "(conditional.x = TRUE) is not scored: lavaan scores only its ",
      "complete cases, and the scores computed here from each case's ",
      "observed items assume conditional.x = FALSE. Refit with ",
      "conditional.x = FALSE to score."), arg = NULL)
  }
  # Without a mean structure lavPredict() takes the means of the observed
  # variables from lav_lisrel_nu0(), which, when the model has no exogenous
  # covariate, sets every intercept to the sample mean, that of an observed
  # variable held as a dummy latent variable included, whose latent
  # variable already carries the mean: the variable's mean is counted
  # twice, which shifts by a constant every score whose weights involve the
  # variable (the regression score of every factor related to it, the
  # Bartlett score of a factor that predicts it). The scores of such a fit
  # are computed here, with the means of .fs_moments() and the formulas
  # lavPredict() uses.
  nu_twice <- !categorical && !conditional_x && any(vapply(est, function(e) {
    is.null(e$nu) && length(e$dummy_y) > 0L && !length(e$dummy_x)
  }, logical(1)))
  # In a model with regressions (a B matrix), lavPredict() gives an observed
  # variable that lavaan holds among the latent variables but that is
  # regressed on nothing (f ~~ z, or a covariate of a factor that also
  # covaries with another variable) a row of zeros in its loading matrix,
  # which drops the variable's covariances with the factors; the scores of
  # such a fit are computed here, with the variable conditioned on
  # (.fs_dummies()). Auxiliary variables are set aside anyway.
  lam_zero <- character(0)
  for (e in est) {
    if (!is.null(e$beta)) {
      lam_zero <- union(lam_zero, setdiff(e$dummy_y, c(e$dummy_r, aux)))
    }
  }
  lam_zero <- if (categorical) character(0) else lam_zero
  if (length(lam_zero) && conditional_x) {
    .dmar_refuse_argument(c(
      "A fit conditional on exogenous covariates (conditional.x = ",
      "TRUE) is not scored when it has regressions and ",
      if (length(lam_zero) == 1L) "an observed variable" else
        "observed variables",
      " that lavaan holds among the latent variables but regresses on ",
      "nothing (", paste(lam_zero, collapse = ", "), "): ",
      "lavaan::lavPredict() then gives such a variable a row of zeros ",
      "among the loadings, which leaves its covariances with the factors ",
      "out of the scores, and the scores computed here assume ",
      "conditional.x = FALSE. Refit with conditional.x = FALSE to ",
      "score."), arg = NULL)
  }
  # The scores are computed here, not by lavPredict(), with auxiliary
  # variables set aside, for a fit that keeps incomplete cases without
  # FIML, for a fit whose means lavaan counts twice, and for one whose
  # loadings lavaan gives a row of zeros.
  own <- length(aux) > 0L || per_case || nu_twice || length(lam_zero) > 0L
  # In a fit of complete rows whose outcome's mean lavaan counts twice, and
  # nothing else, lavaan's Bartlett score of a factor without items of its
  # own (formed through the loadings of the factors it measures) is off by
  # a constant only, which is removed below.
  ho_lavaan <- nu_twice && !per_case && !length(aux) && !length(lam_zero)
  checked <- if (!is.null(data)) .fs_data_check(fit, data) else NULL
  norig <- as.integer(unlist(lavaan::lavInspect(fit, "norig")))
  n <- if (is.null(checked)) max(c(sum(norig), unlist(index))) else checked$n

  bad <- .fs_inadmissible(est, glab, lv, ordered, aux)
  # lavaan computes no EBM or ML score in a group whose solution has a
  # negative residual variance (lav_predict_eta_ebm_ml()), including that of
  # an observed variable regressed on others, which its residual covariance
  # matrix takes from Psi (`theta_eff`).
  neg_theta <- vapply(est, function(e) any(diag(e$theta_eff) < 0),
                      logical(1))
  if (length(bad)) {
    .dmar_warn(c(
      "The solution is inadmissible (", paste(bad, collapse = "; "),
      "). ", if (categorical && any(neg_theta)) {
        paste0("lavaan computes no EBM or ML score for a solution ",
               "with a negative residual variance, so the cases ",
               if (ngroups > 1L) "of that group " else "",
               "have no scores; respecify the model.")
      } else if (categorical) {
        paste0("The EBM scores take the factors' covariance matrix as ",
               "their prior and both scores take the loadings and ",
               "thresholds, so the scores and their standard errors ",
               "inherit the improper solution; respecify the model ",
               "before using them.")
      } else {
        paste0("The regression weights use the model implied ",
               "covariance matrix and the Bartlett weights the ",
               "loadings and the inverse of the residual covariance ",
               "matrix, so the scores, their standard errors, and ",
               "their determinacies inherit the improper solution; ",
               "respecify the model before using them.")
      }), "result_qualified")
  }
  # lavaan's EBM and ML scores of ordered items leave the residual
  # covariances of the observed variables out of each case's likelihood
  # (.fs_residual_pairs()).
  if (categorical) {
    rc <- .fs_residual_pairs(est, lv, glab, skip = neg_theta)
    if (length(rc)) {
      shown <- utils::head(rc, 6L)
      .dmar_warn(c(
        "lavaan's EBM and ML scores of ordered items take each case's ",
        "likelihood as the product of its observed variables' ",
        "probabilities given the factors, which leaves out the ",
        "residual covariance", if (length(rc) > 1L) "s " else " ",
        paste(shown, collapse = ", "),
        if (length(rc) > length(shown)) {
          sprintf(", and %d more", length(rc) - length(shown))
        } else "",
        ": the scores and their standard errors are those of the ",
        "model without ", if (length(rc) > 1L) "them" else "it", "."),
        "result_qualified")
    }
  }
  scored_from <- if (length(aux)) "indicators" else "fit"
  if (length(aux)) {
    several <- length(aux) > 1L
    message("Auxiliary variable", if (several) "s " else " ",
            paste(aux, collapse = ", "), if (several) " enter" else " enters",
            " the fit as saturated correlates. The scores use the fit's ",
            "estimates but are computed from the factor indicators alone, ",
            "since conditioning on an auxiliary that is correlated with ",
            "the item residuals would remove true score variance.")
  }

  items <- .fs_items(fit, lv, ngroups)
  K <- length(lv)
  na_num <- matrix(NA_real_, n, K, dimnames = list(NULL, lv))
  reg <- reg_se <- bart <- bart_se <- det_r <- det_b <- na_num
  i_obs <- i_all <- pat <- matrix(NA_integer_, n, K,
                                  dimnames = list(NULL, lv))
  grp <- rep(NA_character_, n)
  reason <- det_why <- matrix("", n, K, dimnames = list(NULL, lv))
  pattern_tables <- vector("list", ngroups)
  determinacy <- vector("list", ngroups)
  scored <- logical(n)

  why_empty <- "no observed item"
  # A factor none of whose items the case answered: its regression score is
  # the mean of its distribution given whatever else the case has observed
  # (other factors' items, a covariate that predicts it, an outcome it
  # predicts), and the factor's mean when none of that bears on it.
  why_none <- paste0("no item of this factor observed: no ",
                     methods[["bartlett"]], " score, and the ",
                     methods[["regression"]], " score comes only from ",
                     "the case's other observed variables")
  why_none_na <- paste0("no item of this factor observed: no ",
                        methods[["bartlett"]], " score, and no ",
                        methods[["regression"]], " score",
                        if (own) "" else " from lavaan")
  why_extreme <- paste0("every observed item of this factor is at an end ",
                        "of its scale: no finite ML score")
  # The ML score of ordered items (.fs_ml_finite()) and the Bartlett score
  # of continuous items when the observed items do not determine it.
  why_flat <- paste0("the observed items do not separate this factor from ",
                     "another: no ", methods[["bartlett"]], " score")
  why_joint <- paste0("no finite ML score: the likelihood rises without ",
                      "bound as this factor moves with others that share ",
                      "its items")
  why_hi_ml <- paste0("no ML score for a factor without items of its own: ",
                      "the likelihood of the responses does not depend on it")
  why_neg_theta <- paste0("no EBM or ML score: lavaan does not score a ",
                          "solution with a negative residual variance")
  why_det_ord <- "determinacy is not computed for ordered items"
  why_det_cx <- paste0("determinacy is not computed for a fit conditional ",
                       "on exogenous covariates")
  why_det_hi <- paste0("Bartlett determinacy is not computed for a factor ",
                       "without items of its own")
  # With auxiliary variables set aside the scores are computed here from the
  # indicators, and a higher order factor gets no Bartlett score: lavaan
  # forms it through the loadings of the factors it measures, which is not
  # repeated here. Nor does it in the other scores computed here, but for a
  # fit whose only fault is lavaan's double count (`ho_lavaan`).
  why_hi_aux <- paste0("no Bartlett score for a factor without items of ",
                       "its own when the scores are computed from the ",
                       if (length(aux)) {
                         "indicators (auxiliary variables set aside)"
                       } else {
                         paste(fit_desc[2L], "rather than by lavaan")
                       },
                       if (!length(aux) && !per_case && length(lam_zero)) {
                         paste0(", whose loadings give ",
                                paste(lam_zero, collapse = ", "),
                                ", regressed on nothing, a row of zeros")
                       } else "")
  pred_b_ho <- NULL

  if (!own) {
    pred_r <- .fs_lavpredict(fit_lp, methods[["regression"]], lv, model_data)
    pred_b <- .fs_lavpredict(fit_lp, methods[["bartlett"]], lv, model_data)
  }
  # lavaan's ML scores of ordered items with its nlminb optimizer, computed
  # only when needed (see the ordered items below).
  pred_nl <- NULL

  for (g in seq_len(ngroups)) {
    ids <- index[[g]]
    x_g <- model_data[[g]]
    vars <- setdiff(colnames(x_g), aux)
    x_s <- x_g[, vars, drop = FALSE]
    obs <- !is.na(x_s)
    pats <- .fs_patterns(obs)
    grp[ids] <- glab[g]
    scored[ids] <- TRUE
    items_g <- lapply(items[[g]], intersect, vars)
    n_items <- vapply(items_g, length, integer(1))
    io <- vapply(items_g, function(it) {
      as.integer(rowSums(obs[, it, drop = FALSE]))
    }, integer(nrow(obs)))
    io <- matrix(io, nrow(obs), K, dimnames = list(NULL, lv))
    i_obs[ids, ] <- io
    i_all[ids, ] <- matrix(n_items, length(ids), K, byrow = TRUE)
    pat[ids, ] <- pats$id

    # Model matrices restricted to the variables scored from (all of the
    # fit's observed variables, or the indicators when auxiliaries are set
    # aside), for the determinacy and, with auxiliaries, the scores.
    linear <- !categorical && !conditional_x
    if (linear || own) {
      mom <- .fs_moments(est[[g]], x_g)
      keep <- rownames(mom$lambda) %in% vars
      mom$lambda <- mom$lambda[keep, , drop = FALSE]
      mom$theta <- mom$theta[keep, keep, drop = FALSE]
      mom$sigma <- mom$sigma[keep, keep, drop = FALSE]
      mom$ey <- mom$ey[keep]
      lam_vars <- rownames(mom$lambda)
    }
    n_pat <- nrow(pats$observed)
    info <- vector("list", n_pat)
    if (linear) {
      for (p in seq_len(n_pat)) {
        o <- pats$observed[p, match(lam_vars, vars)]
        info[[p]] <- .fs_pattern_info(mom, o)
      }
      full <- .fs_pattern_info(mom, rep(TRUE, length(lam_vars)))
      no_own <- colSums(mom$lambda[, lv, drop = FALSE] != 0) == 0
      det_full_b <- full$det_bart[lv]
      det_full_b[no_own] <- NA_real_
      determinacy[[g]] <- data.frame(
        group = glab[g], factor = lv,
        determinacy_regression = unname(full$det_reg[lv]),
        determinacy_bartlett = unname(det_full_b),
        min_score_correlation = 2 * unname(full$det_reg[lv])^2 - 1,
        stringsAsFactors = FALSE)
    } else {
      no_own <- colSums(est[[g]]$lambda_eff[, lv, drop = FALSE] != 0) == 0
      determinacy[[g]] <- data.frame(
        group = glab[g], factor = lv, determinacy_regression = NA_real_,
        determinacy_bartlett = NA_real_, min_score_correlation = NA_real_,
        stringsAsFactors = FALSE)
    }

    if (own) {
      # Scores from the indicator block (or after two-stage estimation from
      # all of the fit's observed variables), pattern by pattern, with the
      # fit's estimates; the closed forms are those lavPredict() applies.
      yc <- sweep(x_s[, lam_vars, drop = FALSE], 2L, mom$ey[lam_vars])
      r_s <- r_e <- b_s <- b_e <- matrix(NA_real_, nrow(x_s), K,
                                         dimnames = list(NULL, lv))
      for (p in seq_len(n_pat)) {
        rows <- which(pats$id == p)
        ip <- info[[p]]
        if (is.null(ip$reg_w)) next
        o <- pats$observed[p, match(lam_vars, vars)]
        y <- yc[rows, o, drop = FALSE]
        fr <- sweep(y %*% t(ip$reg_w), 2L, mom$eeta, "+")
        r_s[rows, ] <- fr[, lv, drop = FALSE]
        r_e[rows, ] <- matrix(ip$reg_se[lv], length(rows), K, byrow = TRUE)
        if (!is.null(ip$bart_w)) {
          fb <- sweep(y %*% t(ip$bart_w), 2L, mom$eeta, "+")
          b_s[rows, ] <- fb[, lv, drop = FALSE]
          b_e[rows, ] <- matrix(ip$bart_se[lv], length(rows), K,
                                byrow = TRUE)
        }
      }
      reg[ids, ] <- r_s
      reg_se[ids, ] <- r_e
      bart[ids, ] <- b_s
      bart_se[ids, ] <- b_e
      if (ho_lavaan && any(no_own)) {
        # lavaan's Bartlett score of a factor without items of its own, with
        # the constant its double count adds removed: with the observed
        # variables centered at the means of the rows the fit used, as a
        # fit without a mean structure has them, the scores of those
        # (complete) rows average to the factor's mean. Were the mean
        # counted once, the centering would change nothing.
        if (is.null(pred_b_ho)) {
          pred_b_ho <- .fs_lavpredict(fit_lp, "Bartlett", lv, model_data)
        }
        b_ho <- pred_b_ho[[g]]$score[, no_own, drop = FALSE]
        b_ho <- sweep(b_ho, 2L, colMeans(b_ho) - mom$eeta[lv][no_own])
        bart[ids, no_own] <- b_ho
        bart_se[ids, no_own] <- pred_b_ho[[g]]$se[, no_own, drop = FALSE]
      } else {
        for (k in which(no_own)) {
          reason <- .fs_add_reason(reason, ids[rowSums(obs) > 0L], k,
                                   why_hi_aux)
        }
      }
    } else {
      reg[ids, ] <- pred_r[[g]]$score
      reg_se[ids, ] <- pred_r[[g]]$se
      bart[ids, ] <- pred_b[[g]]$score
      bart_se[ids, ] <- pred_b[[g]]$se
    }

    if (linear) {
      det_r[ids, ] <- do.call(rbind, lapply(info, function(ip) {
        ip$det_reg[lv]
      }))[pats$id, , drop = FALSE]
      db <- do.call(rbind, lapply(info, function(ip) {
        ip$det_bart[lv]
      }))[pats$id, , drop = FALSE]
      db[, no_own] <- NA_real_
      det_b[ids, ] <- db
      # Without a Bartlett score (auxiliaries set aside) its reason above
      # already covers the determinacy.
      if (!own || ho_lavaan) det_why[ids, no_own] <- why_det_hi
    } else {
      det_why[ids, ] <- if (categorical) why_det_ord else why_det_cx
    }

    # Continuous items: a pattern whose observed items do not separate a
    # factor from another (an item shared by two factors and nothing else
    # observed of either) has no Bartlett score for it. lavaan's generalized
    # inverse returns the minimum norm solution there, which is not an
    # estimate of the factor's value (.fs_identified(), .fs_pattern_info()).
    if (!categorical) {
      lam_g <- est[[g]]$lambda_eff[, lv, drop = FALSE]
      lam_g <- lam_g[rownames(lam_g) %in% vars, , drop = FALSE]
      at <- match(rownames(lam_g), vars)
      for (p in seq_len(n_pat)) {
        lo <- lam_g[pats$observed[p, at], , drop = FALSE]
        flat <- colSums(lo != 0) > 0 & !.fs_identified(lo)
        if (!any(flat)) next
        hit <- ids[pats$id == p]
        bart[hit, flat] <- bart_se[hit, flat] <- det_b[hit, flat] <- NA_real_
        for (k in which(flat)) {
          reason <- .fs_add_reason(reason, hit, k, why_flat)
        }
      }
    }

    # Ordered items: the ML score of a factor without items of its own, and
    # of a factor whose likelihood has no finite maximum in a case, is not
    # an estimate (.fs_ml_finite()); lavaan returns where its ridge stopped.
    if (categorical && neg_theta[g]) {
      for (k in seq_len(K)) {
        reason <- .fs_add_reason(reason, ids[rowSums(obs) > 0L], k,
                                 why_neg_theta)
      }
    } else if (categorical) {
      for (k in which(no_own)) {
        hit <- ids[rowSums(obs) > 0L]
        bart[hit, k] <- bart_se[hit, k] <- NA_real_
        reason <- .fs_add_reason(reason, hit, k, why_hi_ml)
      }
      ml_code <- .fs_ml_finite(x_s, est[[g]], lv, ordered, !no_own)
      why_ml <- c(extreme = why_extreme, flat = why_flat, joint = why_joint)
      for (k in seq_len(K)) {
        for (cd in names(why_ml)) {
          hit <- ids[ml_code[, k] == cd]
          bart[hit, k] <- bart_se[hit, k] <- NA_real_
          reason <- .fs_add_reason(reason, hit, k, why_ml[[cd]])
        }
      }
      # A case with a factor whose ML score is not finite can stop lavaan's
      # default (BFGS) optimizer on its iteration limit, and lavaan then
      # returns no ML score for any factor of the case. The other factors'
      # ML scores are finite and unique (.fs_ml_finite()): a factor that
      # shares no observed item with the divergent one is separate from it
      # in the likelihood, and one that shares items with it has, in the
      # limit, the maximum of the likelihood of the case's other items,
      # since the shared items' probabilities tend to 1. lavaan's nlminb
      # optimizer reaches either and supplies it.
      divergent <- rowSums(ml_code == "extreme" | ml_code == "joint") > 0L
      b_g <- bart[ids, , drop = FALSE]
      need <- divergent & ml_code == "" &
        matrix(!no_own, length(ids), K, byrow = TRUE) & io > 0L &
        !is.finite(b_g)
      if (any(need)) {
        if (is.null(pred_nl)) {
          pred_nl <- .fs_lavpredict(fit_lp, "ML", lv, model_data,
                                    optim_method = "nlminb")
        }
        e_g <- bart_se[ids, , drop = FALSE]
        b_g[need] <- pred_nl[[g]]$score[need]
        e_g[need] <- pred_nl[[g]]$se[need]
        bart[ids, ] <- b_g
        bart_se[ids, ] <- e_g
      }
    }

    # A factor with none of its items observed has no Bartlett score; the
    # regression score is kept and says where it comes from.
    for (k in seq_len(K)) {
      hit <- ids[io[, k] == 0L & rowSums(obs) > 0L]
      bart[hit, k] <- NA_real_
      bart_se[hit, k] <- NA_real_
      det_b[hit, k] <- NA_real_
      fin <- is.finite(reg[hit, k])
      reason <- .fs_add_reason(reason, hit[fin], k, why_none)
      reason <- .fs_add_reason(reason, hit[!fin], k, why_none_na)
    }
    # A row with no observed item has no score of either kind.
    empty <- ids[rowSums(obs) == 0L]
    reg[empty, ] <- reg_se[empty, ] <- NA_real_
    bart[empty, ] <- bart_se[empty, ] <- NA_real_
    det_r[empty, ] <- det_b[empty, ] <- NA_real_
    reason[empty, ] <- why_empty
    det_why[empty, ] <- ""

    # The per-pattern table.
    miss_names <- apply(pats$observed, 1L, function(o) {
      paste(vars[!o], collapse = ", ")
    })
    pattern_tables[[g]] <- do.call(rbind, lapply(seq_len(n_pat), function(p) {
      o <- pats$observed[p, ]
      ip <- info[[p]]
      d_r <- if (linear) ip$det_reg[lv] else rep(NA_real_, K)
      d_b <- if (linear) ip$det_bart[lv] else rep(NA_real_, K)
      s_r <- if (linear) ip$reg_se[lv] else rep(NA_real_, K)
      s_b <- if (linear) ip$bart_se[lv] else rep(NA_real_, K)
      io_p <- vapply(items_g, function(it) sum(o[match(it, vars)]),
                     numeric(1))
      # lavaan's Bartlett score of a factor without items of its own (formed
      # through the loadings of the factors it measures) has one standard
      # error per pattern, which the rows carry.
      if (linear && (!own || ho_lavaan) && any(no_own)) {
        rows_p <- ids[pats$id == p]
        s_b[no_own] <- vapply(which(no_own), function(k) {
          v <- bart_se[rows_p, k]
          v <- v[is.finite(v)]
          if (length(v)) v[1L] else NA_real_
        }, numeric(1))
      }
      d_b[no_own | io_p == 0] <- NA_real_
      s_b[io_p == 0] <- NA_real_
      if (!any(o)) d_r <- s_r <- rep(NA_real_, K)
      data.frame(group = glab[g], pattern = p, factor = lv,
                 n = sum(pats$id == p),
                 items_observed = as.integer(io_p),
                 items = unname(n_items),
                 missing_variables = miss_names[p],
                 determinacy_regression = unname(d_r),
                 determinacy_bartlett = unname(d_b),
                 min_score_correlation = 2 * unname(d_r)^2 - 1,
                 regression_se = unname(s_r),
                 bartlett_se = unname(s_b),
                 stringsAsFactors = FALSE)
    }))
  }

  # A score without a standard error says so; for continuous items that
  # means a score variance below 1e-5, which lavPredict() reports as NA and
  # the scores computed here (from the indicators, or after two-stage
  # estimation) treat the same way. Only what lavaan computed is attributed
  # to lavaan.
  floor_note <- if (categorical) "" else " (a score variance below 1e-5)"
  source_note <- if (own) "" else " from lavaan"
  why_missing <- if (own) {
    "score not computable from the observed indicators"
  } else {
    "score unavailable from lavaan"
  }
  for (k in seq_len(K)) {
    hit <- which(is.finite(reg[, k]) & !is.finite(reg_se[, k]))
    reason <- .fs_add_reason(reason, hit, k, paste0(
      "no ", methods[["regression"]], " standard error", source_note,
      floor_note))
    hit <- which(is.finite(bart[, k]) & !is.finite(bart_se[, k]))
    reason <- .fs_add_reason(reason, hit, k, paste0(
      "no ", methods[["bartlett"]], " standard error", source_note,
      floor_note))
    hit <- which(scored & !nzchar(reason[, k]) &
                   (!is.finite(reg[, k]) | !is.finite(bart[, k])))
    reason <- .fs_add_reason(reason, hit, k, why_missing)
  }
  # Determinacy that is not computed for a whole factor or fit comes last.
  for (txt in setdiff(unique(as.vector(det_why)), "")) {
    for (k in seq_len(K)) {
      reason <- .fs_add_reason(reason, which(det_why[, k] == txt), k, txt)
    }
  }

  # Rows the fit did not score: deleted listwise, deleted for a missing
  # exogenous covariate, without a group, or in a group the fit left out.
  # With `data` the reason is identified and the items answered counted.
  unscored <- which(!scored)
  if (length(unscored)) {
    miss <- opts$missing
    # One reason each, without the "; " that separates reasons. lavaan has
    # no FIML for ordered items; pairwise deletion keeps the case, and the
    # EBM and ML scores use its observed items.
    why_lw <- if (categorical) {
      paste0("deleted listwise by the fit (with ordered items, a fit with ",
             "missing = \"pairwise\" scores it)")
    } else {
      paste0("deleted listwise by the fit (a fit by FIML, missing = ",
             "\"ml\", scores it)")
    }
    # A row missing an exogenous covariate is deleted by any fit with
    # fixed.x = TRUE (lavaan's default) but one with missing = "ml.x", and
    # by any fit conditional on the covariates, whatever its treatment of
    # missing values: FIML with missing = "ml" deletes it too, and so does
    # pairwise deletion of ordered items. Only missing = "ml.x" (refused
    # with conditional.x = TRUE) or fixed.x = FALSE keeps it, for
    # continuous items; with ordered items lavaan keeps it only in a fit
    # not conditional on the covariates, which lavPredict() cannot score
    # (.fs_check_dummies()).
    deletes_x <- length(ov_x) > 0L && !identical(miss, "ml.x") &&
      (isTRUE(opts$fixed.x) || conditional_x || categorical)
    why_cov <- if (categorical) {
      paste0("deleted by the fit for a missing exogenous covariate (with ",
             "ordered items, no fit that lavaan can score keeps it)")
    } else if (conditional_x) {
      paste0("deleted by the fit for a missing exogenous covariate ",
             "(conditional.x = TRUE), which missing = \"ml.x\" with ",
             "conditional.x = FALSE scores")
    } else {
      paste0("deleted by the fit for a missing exogenous covariate ",
             "(fixed.x = TRUE), which missing = \"ml.x\" scores")
    }
    if (!is.null(checked)) {
      d <- checked$data
      for (r in unscored) {
        g_r <- if (!is.null(checked$group)) {
          as.character(d[[checked$group]][r])
        } else {
          NA_character_
        }
        grp[r] <- g_r
        g_idx <- if (ngroups > 1L) match(g_r, glab) else 1L
        it_r <- if (!is.na(g_idx)) items[[g_idx]] else
          stats::setNames(lapply(lv, function(f) {
            unique(unlist(lapply(items, `[[`, f)))
          }), lv)
        i_all[r, ] <- vapply(it_r, length, integer(1))
        i_obs[r, ] <- vapply(it_r, function(it) {
          as.integer(sum(!is.na(unlist(d[r, it]))))
        }, integer(1))
        why <- if (r %in% checked$unassigned) {
          paste0("not used by the fit: the group variable ",
                 checked$group, " is missing")
        } else if (r %in% checked$unselected) {
          "not used by the fit: a group the fit leaves out (group.label)"
        } else if (deletes_x && anyNA(unlist(d[r, ov_x]))) {
          why_cov
        } else if (anyNA(unlist(d[r, ov]))) {
          why_lw
        } else {
          "not used by the fit"
        }
        reason[r, ] <- why
      }
    } else {
      # Without `data` the row's values are unknown, so a row deleted for a
      # missing exogenous covariate cannot be told from one deleted for a
      # missing item. When the fit deletes both kinds, the reason names
      # what scores each; telling it that missing = "ml" (or, with ordered
      # items, "pairwise") scores it would be false for the first kind.
      x_remedy <- if (conditional_x) {
        "missing = \"ml.x\" with conditional.x = FALSE"
      } else {
        "missing = \"ml.x\""
      }
      why <- if (identical(miss, "listwise") && ngroups == 1L) {
        if (!deletes_x) {
          why_lw
        } else if (categorical) {
          paste0("deleted listwise by the fit (with ordered items, a fit ",
                 "with missing = \"pairwise\" scores a row with a missing ",
                 "item, no fit that lavaan can score keeps one with a ",
                 "missing exogenous covariate, and supplying data ",
                 "identifies the reason)")
        } else {
          paste0("deleted listwise by the fit (a fit by FIML, missing = ",
                 "\"ml\", scores a row with a missing item, ",
                 x_remedy, " one with a missing ",
                 "exogenous covariate, and supplying data identifies the ",
                 "reason)")
        }
      } else if (deletes_x) {
        paste0("not used by the fit, which deletes a row with a missing ",
               "exogenous covariate (", if (categorical) {
                 "with ordered items no fit that lavaan can score keeps "
               } else {
                 paste0(x_remedy, " scores ")
               }, "such a row, and supplying data identifies the reason)")
      } else {
        "not used by the fit"
      }
      reason[unscored, ] <- why
      if (ngroups == 1L) {
        i_all[unscored, ] <- matrix(vapply(items[[1L]], length, integer(1)),
                                    length(unscored), K, byrow = TRUE)
      }
    }
  }

  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  for (g in which(!is.na(bartlett_issues))) {
    ids <- index[[g]]
    bart[ids, ] <- bart_se[ids, ] <- det_b[ids, ] <- NA_real_
    determinacy[[g]]$determinacy_bartlett <- NA_real_
    pattern_tables[[g]]$determinacy_bartlett <- NA_real_
    pattern_tables[[g]]$bartlett_se <- NA_real_
    for (k in seq_len(K)) {
      reason <- .fs_add_reason(reason, ids, k, bartlett_issues[g])
    }
  }
  # The score columns carry the prefix `.derived_`, so that a score joined
  # to the observed data cannot be taken for an observed variable or
  # collide with one. check.names = FALSE keeps the leading dot.
  out <- data.frame(
    row = rep(seq_len(n), times = K),
    group = rep(grp, times = K),
    factor = rep(lv, each = n),
    .derived_regression = as.vector(reg),
    .derived_regression_se = as.vector(reg_se),
    .derived_regression_ci_lower = as.vector(reg - z * reg_se),
    .derived_regression_ci_upper = as.vector(reg + z * reg_se),
    .derived_bartlett = as.vector(bart),
    .derived_bartlett_se = as.vector(bart_se),
    .derived_bartlett_ci_lower = as.vector(bart - z * bart_se),
    .derived_bartlett_ci_upper = as.vector(bart + z * bart_se),
    items_observed = as.vector(i_obs),
    items = as.vector(i_all),
    pattern = as.vector(pat),
    determinacy_regression = as.vector(det_r),
    determinacy_bartlett = as.vector(det_b),
    reason = as.vector(reason),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  # The eight score columns, named by what each holds, so that code can
  # look a column up by its role (attr(x, "derived_columns")[["bartlett"]]).
  score_columns <- c("regression", "regression_se", "regression_ci_lower",
                     "regression_ci_upper", "bartlett", "bartlett_se",
                     "bartlett_ci_lower", "bartlett_ci_upper")
  derived_columns <- stats::setNames(paste0(".derived_", score_columns),
                                     score_columns)
  class(out) <- c("dmar_factor_scores", "data.frame")
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  patterns <- do.call(rbind, pattern_tables)
  rownames(patterns) <- NULL
  determinacy <- do.call(rbind, determinacy)
  rownames(determinacy) <- NULL
  attr(out, "methods") <- methods
  attr(out, "derived_columns") <- derived_columns
  attr(out, "bartlett_issues") <- bartlett_issues
  attr(out, "missing") <- opts$missing
  attr(out, "estimator") <- opts$estimator.orig %||% opts$estimator
  attr(out, "patterns") <- patterns
  attr(out, "determinacy") <- determinacy
  attr(out, "scored_from") <- scored_from
  attr(out, "auxiliary") <- aux
  attr(out, "ngroups") <- ngroups
  attr(out, "lavaan_version") <- as.character(utils::packageVersion("lavaan"))
  out
}

#' @rdname factor_scores
#' @param x A \code{dmar_factor_scores} table returned by
#'   \code{factor_scores()}.
#' @param n For \code{print()}, the number of rows shown for each factor
#'   (default 3); the pattern table shows up to \code{4 * n} rows.
#'   \code{\link{knit_print.dmar_tbl}} and \code{\link{as_kable}} show
#'   the same summary with the default.
#' @param digits For \code{print()}, the decimal places of the scores,
#'   their standard errors and limits, and the determinacies (default
#'   \code{getOption("dmar.digits", 3)}; see \code{\link{dmar_tbl}}).
#' @param \dots Unused.
#' @export
print.dmar_factor_scores <- function(x, n = 3L,
                                     digits = getOption("dmar.digits", 3L),
                                     ...) {
  # A subset of the table (x[rows, cols]) keeps the class but not the
  # attributes the summary needs, so it prints as an ordinary DMAR table.
  if (!.factor_scores_summary_ok(x)) return(invisible(NextMethod()))
  .dmar_tbl_write(x, digits = digits, print_args = list(right = TRUE),
                  view_args = list(n = n, digits = digits))
}

# Every display of a factor_scores() table (the console, knit_print(),
# as_kable(); see .dmar_tbl_display() and .dmar_tbl_footer() in
# R/dmar_output_helpers.R) is a summary rather than the whole table: two
# lines that say what was scored and from which fit, the first `n` rows
# of each factor with their intervals, the reasons values are missing or
# borrowed, the determinacy of each missing data pattern, and how the
# intervals are to be read. A subset that has lost the attributes the
# summary reads displays as an ordinary DMAR table.
.factor_scores_summary_ok <- function(x) {
  !is.null(attr(x, "methods")) &&
    all(c("factor", "row", ".derived_regression", ".derived_bartlett",
          "reason") %in% names(x))
}

#' @exportS3Method .dmar_tbl_display dmar_factor_scores
#' @noRd
.dmar_tbl_display.dmar_factor_scores <- function(
    x, n = 3L, digits = getOption("dmar.digits", 3L), ...) {
  if (!.factor_scores_summary_ok(x)) return(NextMethod())
  methods <- attr(x, "methods")
  lv <- unique(x$factor)
  fmt <- function(v) {
    out <- .format_dmar_value(v, digits, whole = FALSE)
    out[is.na(out)] <- "NA"
    out
  }
  ci <- function(lo, hi) {
    ifelse(is.na(lo) | is.na(hi), "NA",
           sprintf("[%s, %s]", fmt(lo), fmt(hi)))
  }
  shown <- unlist(lapply(lv, function(f) {
    utils::head(which(x$factor == f), n)
  }))
  h <- as.data.frame(x)[shown, , drop = FALSE]
  level <- sprintf("%s%% CI", .format_conf_pct(attr(x, "conf_level")))
  disp <- data.frame(row = h$row, stringsAsFactors = FALSE)
  if (any(!is.na(x$group))) disp$group <- ifelse(is.na(h$group), "NA", h$group)
  disp$factor <- h$factor
  disp[[methods[["regression"]]]] <- fmt(h$.derived_regression)
  disp[["se"]] <- fmt(h$.derived_regression_se)
  disp[[level]] <- ci(h$.derived_regression_ci_lower,
                      h$.derived_regression_ci_upper)
  disp[[methods[["bartlett"]]]] <- fmt(h$.derived_bartlett)
  disp[["se "]] <- fmt(h$.derived_bartlett_se)
  disp[[paste0(level, " ")]] <- ci(h$.derived_bartlett_ci_lower,
                                   h$.derived_bartlett_ci_upper)
  .as_dmar_tbl(disp)
}

#' @exportS3Method .dmar_tbl_footer dmar_factor_scores
#' @noRd
.dmar_tbl_footer.dmar_factor_scores <- function(
    x, n = 3L, digits = getOption("dmar.digits", 3L), ...) {
  if (!.factor_scores_summary_ok(x)) return(NextMethod())
  lines <- .factor_scores_footer(x, n, digits)
  # The confidence level (from the default method) closes the summary.
  c(lines, NextMethod())
}

# The lines of a factor_scores() summary, named by block (see
# .dmar_tbl_footer()): the two "header" lines that introduce the table, the
# "table_note" under it, the counts of each reason a value is missing or
# borrowed, the determinacy of each missing data pattern (its rows as a
# small "_table"), and the reading of the intervals.
.factor_scores_footer <- function(x, n = 3L,
                                  digits = getOption("dmar.digits", 3L)) {
  methods <- attr(x, "methods")
  lv <- unique(x$factor)
  ngroups <- attr(x, "ngroups") %||% 1L
  aux <- attr(x, "auxiliary")
  grouped <- any(!is.na(x$group))
  rows <- length(unique(x$row))
  fmt <- function(v) {
    out <- .format_dmar_value(v, digits, whole = FALSE)
    out[is.na(out)] <- "NA"
    out
  }
  lines <- .dmar_footer_add(character(0), c(
    sprintf("Factor scores: %s and %s scores, %d %s (%s), %d rows",
            methods[["regression"]], methods[["bartlett"]], length(lv),
            if (length(lv) == 1L) "factor" else "factors",
            paste(lv, collapse = ", "), rows),
    sprintf("Fit: %s, missing = \"%s\"%s; scored from %s",
            attr(x, "estimator"), attr(x, "missing"),
            if (ngroups > 1L) sprintf(", %d groups", ngroups) else "",
            if (identical(attr(x, "scored_from"), "indicators")) {
              sprintf("the indicators (auxiliary %s set aside)",
                      paste(aux, collapse = ", "))
            } else {
              "the fit"
            })), "header")
  lines <- .dmar_footer_add(lines, if (isTRUE(n == 1)) {
    "(the first row of each factor)"
  } else {
    sprintf("(the first %d rows of each factor)", n)
  }, "table_note")

  why <- unlist(strsplit(x$reason[nzchar(x$reason)], "; ", fixed = TRUE))
  if (length(why)) {
    tab <- sort(table(why), decreasing = TRUE)
    lines <- .dmar_footer_add(lines, c(
      paste("Missing or borrowed values (row and factor pairs; see the",
            "reason column):"),
      sprintf("  %d: %s", as.integer(tab), names(tab))), "missing_values")
  }

  pt <- attr(x, "patterns")
  if (!is.null(pt) && any(is.finite(pt$determinacy_regression))) {
    show <- utils::head(pt, 4L * n)
    dt <- data.frame(pattern = show$pattern, n = show$n,
                     stringsAsFactors = FALSE)
    if (grouped) dt <- cbind(group = show$group, dt)
    dt$factor <- show$factor
    dt$items <- paste0(show$items_observed, "/", show$items)
    dt[[methods[["regression"]]]] <- fmt(show$determinacy_regression)
    dt[[methods[["bartlett"]]]] <- fmt(show$determinacy_bartlett)
    if (!is.null(show$min_score_correlation))
      dt$min_r <- fmt(show$min_score_correlation)
    lines <- .dmar_footer_add(lines, paste(
      "Determinacy, the correlation of each score with its factor,",
      "by missing data pattern; min_r, the lowest correlation two sets of",
      "scores equally consistent with the model can have (2 x the",
      "regression score's determinacy squared, minus 1; Guttman, 1955):"),
      "determinacy")
    lines <- .dmar_footer_add(
      lines, utils::capture.output(print(dt, row.names = FALSE, right = TRUE)),
      "determinacy_table")
    if (nrow(pt) > nrow(show)) {
      lines <- .dmar_footer_add(lines, sprintf(
        "(%d of %d rows; all of them in attr(x, \"patterns\"))",
        nrow(show), nrow(pt)), "determinacy")
    }
  }

  # One paragraph, broken for the console within 80 columns; a knitted
  # report reads the breaks as spaces.
  notes <- if (identical(methods[["regression"]], "EBM")) {
    c("Intervals treat the fitted parameters as known: EBM intervals are",
      "approximate posterior intervals around the posterior mode, ML intervals",
      "Wald intervals for each case's factor value. For relations among",
      "factors, one-step SEM is preferred to any regression among scores.")
  } else {
    c("Intervals treat the fitted parameters as known: regression intervals",
      "are posterior (empirical Bayes) intervals, Bartlett intervals Wald",
      "intervals for each case's factor value. In a later regression, use",
      "regression scores for a predictor and Bartlett scores for an outcome,",
      "each from its own measurement model (Skrondal & Laake, 2001); for",
      "relations among factors, one-step SEM or structural after measurement",
      "(lavaan::sam()) is preferred to any regression among scores.")
  }
  .dmar_footer_add(lines, paste(notes, collapse = "\n"), "notes")
}

#' @rdname factor_scores
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_factor_scores <- function(x, ...) {
  out <- x
  for (a in setdiff(names(attributes(out)), c("names", "row.names"))) {
    attr(out, a) <- NULL
  }
  class(out) <- "data.frame"
  rownames(out) <- NULL
  out
}

#' @rdname factor_scores
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_factor_scores <- function(x, ...) {
  methods <- attr(x, "methods") %||%
    c(regression = NA_character_, bartlett = NA_character_)
  pt <- attr(x, "patterns")
  cols <- intersect(c(".derived_regression", ".derived_bartlett"), names(x))
  scored <- if (length(cols) && !is.null(x$row)) {
    fin <- Reduce(`|`, lapply(cols, function(v) is.finite(x[[v]])))
    sum(tapply(fin, x$row, any))
  } else {
    NA_integer_
  }
  data.frame(
    n_rows = length(unique(x$row)),
    n_scored = scored,
    n_factors = length(unique(x$factor)),
    n_groups = attr(x, "ngroups") %||% 1L,
    n_patterns = if (is.null(pt)) NA_integer_ else
      nrow(unique(pt[, c("group", "pattern")])),
    method_regression = methods[["regression"]],
    method_bartlett = methods[["bartlett"]],
    estimator = attr(x, "estimator") %||% NA_character_,
    missing = attr(x, "missing") %||% NA_character_,
    scored_from = attr(x, "scored_from") %||% NA_character_,
    conf_level = attr(x, "conf_level") %||% NA_real_,
    stringsAsFactors = FALSE
  )
}
