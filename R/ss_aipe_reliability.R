#' Sample Size Planning for Accuracy in Parameter Estimation for Reliability Coefficients
#'
#' @description
#' Computes the necessary sample size for the confidence interval on a
#' population reliability coefficient (coefficient alpha or coefficient
#' omega, depending on the assumed measurement model) to have expected
#' width no larger than \code{width}, or (when \code{assurance} is
#' supplied) to be no wider than \code{width} with the specified
#' probability. This is the accuracy in parameter estimation (AIPE)
#' counterpart to power based planning for reliability and is the
#' companion of \code{\link{reliability}}. Without an assurance the
#' sample size comes from a closed form for the expected width of the
#' interval \code{type} names: the normal theory variance of coefficient
#' alpha of van Zyl, Neudecker, and Nel (2000) for the parallel model and
#' for the normal theory interval of the tau equivalent model (Terry &
#' Kelley, 2012, Eq. 19), and the delta method variance of the model
#' implied coefficient for the factor analytic interval of the tau
#' equivalent model (equal loadings) and for the congeneric model (coefficient
#' omega). When \code{assurance} is supplied, a
#' Monte Carlo search starting from the closed form finds the smallest
#' \emph{N} delivering the requested assurance (Terry & Kelley, 2012).
#'
#' Coefficient alpha (Guttman, 1945; subsequently popularized by
#' Cronbach, 1951) is planned for under the parallel and tau equivalent
#' (\emph{i.e.}, the so called True Score) measurement models;
#' coefficient omega (McDonald, 1999) is planned for under the congeneric
#' model. The coefficient is the one \code{\link{reliability_alpha}} and
#' \code{\link{reliability_omega}} estimate from the raw item scores: the
#' population, however it is given, is used as a covariance matrix, so
#' items of unequal variance are planned for as they are, and a correlation
#' matrix plans for standardized items.
#'
#' When the coefficient will be estimated from incomplete data by full
#' information maximum likelihood, the interval at a given \emph{N} is
#' wider than it would be with complete data. The \code{fmi} argument
#' plans for that: it is the fraction of missing information anticipated
#' for the coefficient, and the plan grows to carry the information that
#' complete data would have carried (see the section on planning under
#' missing data).
#'
#' @details
#' \strong{The closed forms.} Each closed form takes the sampling variance
#' of the coefficient to be \eqn{V / (N - 1)} and solves
#' \eqn{2 z \sqrt{V / (N - 1)} = w}{2 z sqrt(V / (N - 1)) = w} for
#' \emph{N}, so that \eqn{N = 4 z^2 V / w^2 + 1}, where \emph{w} is
#' \code{width} and \emph{z} the standard normal quantile for
#' \code{conf_level}. For the parallel model \emph{V} is the compound
#' symmetry form \eqn{2 (1 - \alpha)^2 J / (J - 1)}{2 (1 - alpha)^2 J /
#' (J - 1)} of van Zyl et al. (2000), with \emph{J} items, which at a
#' population of parallel items is also the variance of the factor
#' analytic interval, so the parallel model has one closed form for both
#' types. For the normal theory interval of the tau equivalent model
#' \emph{V} is the general variance of alpha of van Zyl et al. (Terry &
#' Kelley, 2012, Eqs. 7, 8, and 19), evaluated on the population covariance
#' matrix as the formula defines it. For the congeneric model
#' \emph{V} is the delta method variance of the model implied omega,
#' \eqn{\omega = (\sum \lambda_j)^2 / ((\sum \lambda_j)^2 + \sum
#' \psi_j^2)}{omega = (sum lambda)^2 / ((sum lambda)^2 + sum psi^2)}, with
#' \eqn{\lambda_j} the loadings (\code{true_factor_loading}) and
#' \eqn{\psi_j^2} the error variances (\code{true_residual_var}): the
#' gradient of omega taken through the inverse of the expected information
#' of the single-factor model at the population (Raykov, 2002; Terry &
#' Kelley, 2012, Eqs. 13 to 16), the variance behind
#' \code{reliability_omega(denominator = "model_implied", ci_method =
#' "ml")}, whose squared standard error is \eqn{V / N}; taking it as
#' \eqn{V / (N - 1)}, as for alpha, adds one case and keeps the planning
#' under missing data below the same for every model. When the population
#' is given by \code{true_cov_matrix} or \code{data}, the single-factor
#' model is first fit to it by maximum likelihood. Terry and Kelley (2012,
#' p. 390) had no closed form for omega: they used the alpha formula as
#' if true score
#' equivalence held only as a starting value, and found the size whose
#' mean interval width met the target by an a priori Monte Carlo
#' simulation. That starting value, taken as the plan for omega, is too
#' large when the loadings differ among standardized items and can be far
#' too small when the item variances differ. The delta method form gives
#' the expected width to first order. The mean width of the factor
#' analytic interval at the planned \emph{N} meets the target where the
#' samples at that \emph{N}
#' stay away from an improper solution and from empirical
#' underidentification, and exceeds it where some do not, whose intervals
#' are far wider, while the median width stays near the target (within 4
#' percent of it in the simulations below, above it in some). Which samples
#' come near depends on the loadings, the number of items, and \emph{N}
#' together, so no bound on the loadings alone marks the populations the
#' closed form serves. With three items the single-factor model is just
#' identified, each squared loading is a ratio of the covariances
#' (\eqn{\lambda_1^2 = \sigma_{12} \sigma_{13} / \sigma_{23}}{lambda_1^2 =
#' sigma_12 sigma_13 / sigma_23}), and a sample covariance near zero makes
#' the solution run away: three loadings of .3 planned for a width of .20
#' take 687 cases, at which the median width is .202 and the mean .264,
#' and 3 percent of samples give no interval (1,000 samples), while at a
#' width of .10 (2743 cases) the mean is on target; loadings of .5, .4,
#' and .3 at .20 (462 cases) average .242, with a median of .206, and
#' three loadings of .4 (468 cases) .211, with a median of .199 (500
#' samples each). With four loadings of .3 the mean
#' width was on target at .20 (527 cases) but averaged .323 at .30 (235
#' cases), where a tenth of the samples give no interval. Where the factor
#' rests on weak loadings the same happens with more items: loadings of
#' .7, .7, .1, and .1 (omega .46) planned for a width of .10 take 974
#' cases, at which the median width is .104 and the mean .181, and four
#' loadings of .15 take 3437 cases, with a median of .098 and a mean of
#' .148. In the other populations simulated for this planner, with four
#' or more items whose loadings share a sign, the mean was within a few
#' percent of the target (for loadings of opposite sign, see below). For
#' a population like those above, plan with an assurance above .5, whose
#' search counts the widths themselves.
#'
#' For the factor analytic interval of the tau equivalent model \emph{V}
#' is the same delta method variance with the loadings held equal, the
#' variance behind \code{reliability_alpha(estimator = "model_implied",
#' ci_method = "ml")}: the single loading's information is the sum of the
#' loading block of the information above, and its gradient the sum of the
#' loadings' gradients. That coefficient is the maximum likelihood
#' estimator under the model, so its variance is no larger than the normal
#' theory variance of alpha, and it is smaller whenever the error
#' variances differ: with loadings .7 and error variances .1, .2, 1, and 2,
#' a width of .10 takes 242 cases with the factor analytic interval and
#' 367 with the normal theory one. With equal error variances the two
#' coincide.
#'
#' \strong{Populations that are refused.} The planning population has to
#' be one the model describes. A population given as
#' \code{true_cov_matrix} must be a covariance matrix (no negative
#' eigenvalue), and when the single-factor model is fit to it, the fit must
#' be proper. Where that model is fit, for the congeneric model and for the
#' tau equivalent model with the factor analytic interval, every error
#' variance in \code{true_residual_var} must be above zero, as the parallel
#' model requires: a negative one is not a variance, and a zero
#' one makes an item free of error, the boundary at which the single-factor
#' solution is improper. The same population given as
#' \code{true_cov_matrix} is refused through its fit. With the normal
#' theory interval of the tau equivalent model nothing is fit: the loading
#' and the error variances only build the covariance matrix the plan is
#' made from, which must be a covariance matrix, as \code{true_cov_matrix}
#' must, and the plan is the one the same matrix given as
#' \code{true_cov_matrix} gives. The congeneric model needs at least three
#' items with a nonzero loading, without which it is not identified at
#' the population
#' and omega has no sampling variance. With the factor analytic interval of
#' the parallel and tau equivalent models, items to which the
#' single-factor model with equal loadings gives a loading of zero share
#' no common factor in that model, and the coefficient that interval
#' estimates is zero, on the boundary of its range, where the closed form
#' does not describe the interval. For a population given by its loading
#' (\code{true_item_covariance} or \code{true_factor_loading}) that is a
#' loading of zero, and for one given as \code{true_cov_matrix} or
#' \code{data} the loading of the maximum likelihood fit. For items of
#' equal variance it happens exactly when their covariances average zero
#' or below (alpha of the population matrix at or below zero); with
#' unequal variances the fit weighs each covariance by
#' the inverse of its two items' variances, and the two can disagree either
#' way (three items of variance 1 that correlate .3 and a fourth of
#' variance 4 that correlates -.25 with them have alpha -0.276 and a
#' fitted coefficient of .315, which the tau equivalent model plans for).
#' The parallel model's closed form is a function of alpha of the
#' population matrix, and it is refused when alpha is at or below zero; a
#' population whose fitted coefficient is nonetheless positive is far from
#' parallel, and the message names the tau equivalent model.
#' A supplied population under \code{model = "Parallel"} must have equal
#' diagonal entries and equal off-diagonal entries, each within
#' \eqn{\sqrt{\epsilon}\max |S|}, where \eqn{\epsilon} is machine
#' precision. This check uses the original item units. A sample covariance
#' will generally depart from an exact parallel population; specify the
#' intended population with \code{true_factor_loading} and
#' \code{true_residual_var} when parallelism is a planning assumption rather
#' than a property of \code{true_cov_matrix}.
#'
#' \strong{Where the congeneric closed form is refused.} The delta method
#' takes omega, a square of the sum of the loadings \eqn{L} over the
#' variance of the composite, to be normal with its first-order variance,
#' which needs \eqn{L} to be large relative to its standard error. When it
#' is not, as when loadings of opposite sign nearly cancel in the
#' unweighted sum, omega behaves near zero as the square of a normal
#' variable, and its first-order variance understates the width of its
#' interval; where the loadings cancel exactly, the gradient of omega
#' vanishes and the closed form gives one case. The closed form is
#' therefore checked at the \emph{N} it gives, or at \code{n_items + 1} cases
#' when it gives fewer, the fewest at which the factor analytic interval can
#' be computed (with fewer, the sample covariance matrix is singular). Write
#' \eqn{\sigma_L}{sigma_L}
#' for the standard error of \eqn{\hat L}{L_hat} there, \eqn{r = \sigma_L
#' / |L|}{r = sigma_L / |L|}, and \eqn{P} for the sum of the error
#' variances. Near zero
#' \eqn{\hat\omega}{omega_hat} is about \eqn{\hat L^2 / P}{L_hat^2 / P}
#' and its delta method standard error about \eqn{2 |\hat L| \sigma_L /
#' P}{2 |L_hat| sigma_L / P}, so with \eqn{\hat L}{L_hat} normal, the
#' interval the estimation functions report, clamped at zero, averages
#' \deqn{\frac{E[\min(4 z |X|, X^2 + 2 z |X|)]}{4 z / r}, \qquad X \sim
#' N(1 / r, 1),}{E[min(4 z |X|, X^2 + 2 z |X|)] / (4 z / r), X ~ N(1 / r,
#' 1),}
#' times the width planned for. That ratio stays below one, the clamped
#' interval making the plan conservative, until \eqn{r} is about 1.5 at 95
#' percent confidence (1.36 at 90 percent and 1.66 at 99), and it exceeds
#' one beyond; the closed form is refused once it does. In simulations of
#' the factor analytic interval at the planned \emph{N} for loadings .7,
#' -.7, .5, and a fourth loading of -.2, -.32, -.35, or -.37 (\eqn{r} of
#' .41, 1.07, 1.52, and 2.00 there, the last two refused), the mean width was
#' .082, .082, .089, and .121 for a target of .10. With a fourth loading of
#' -.4 (omega .0038, \eqn{r = 3.35}) the closed form's 49 cases give a mean
#' width of .22, where the interval needs about 95, and with -.5 the
#' loadings cancel and the delta method gives one case. An assurance plan
#' is not refused: its search counts the widths of the intervals
#' themselves, starting at no fewer than \code{n_items + 1} cases.
#' A congeneric closed form below \code{n_items + 1} is floored at that size
#' before any missing-information inflation, with a warning. This is only
#' a feasibility minimum. It does not establish that the expected width
#' target is met at such a small sample size; an assurance plan is
#' recommended.
#'
#' Whenever the loadings of a congeneric population differ in sign, the
#' planner warns. An item that loads opposite to the others lowers omega of
#' the unweighted sum, as an item keyed in the opposite direction does when
#' it is not reverse scored. The plan is for omega of the sum of the items
#' as they are given; reverse score such items before planning, as before
#' the coefficient is estimated.
#'
#' The parallel and tau equivalent models need no such check. Their one
#' loading enters the covariance matrix only through its square, the
#' common covariance, of which their coefficient is a smooth function, so
#' its variance does not vanish near zero, and there the clamped interval
#' only makes the plan conservative: four parallel items that covary .03
#' (alpha .11) planned for a width of .30 take 362 cases, at which the
#' factor analytic interval averages .236 wide.
#'
#' \strong{The assurance search.} The Monte Carlo search draws samples
#' from the population covariance matrix and, at each candidate sample
#' size, computes the realized confidence interval width with the same
#' machinery the estimation side of the package uses. For
#' \code{type = "Factor Analytic"} the single-factor model is fit by
#' maximum likelihood (equal loadings for the parallel and tau
#' equivalent models, free loadings for the congeneric model) and the
#' interval is the delta method Wald interval on the model implied
#' reliability, the interval of
#' \code{\link{reliability_omega}(denominator = "model_implied",
#' ci_method = "ml")} and of
#' \code{\link{reliability_alpha}(estimator = "model_implied",
#' ci_method = "ml")}. For \code{type = "Normal Theory"} the interval
#' uses the van Zyl, Neudecker, and Nel (2000) closed form standard
#' error for the tau equivalent model and its compound symmetry
#' simplification for the parallel model, the same closed form behind
#' \code{reliability_alpha(ci_method = "ml")}. The congeneric model has
#' no normal theory form, so \code{type = "Normal Theory"} is an error
#' there.
#'
#' The search starts from the closed form size, and with the factor
#' analytic interval at no fewer than \code{n_items + 1} cases, the fewest whose
#' sample covariance matrix the single-factor model can be fit to (a
#' highly reliable population planned for a wide interval can have a
#' closed form below that). When that size falls
#' short of the assurance, the size is raised one case at a time until a
#' pass of \code{G_initial} samples meets it, and a final pass of
#' \code{G} samples confirms the size, again raising it if needed
#' (Terry & Kelley, 2012). When the closed form size already meets the
#' assurance, which happens when the assurance is near one half or below,
#' or when the interval of \code{type} is narrower than the one the closed
#' form describes (as the factor analytic interval of the parallel model
#' can be when the population given as \code{true_cov_matrix} or
#' \code{data} is not parallel), the search moves down: to
#' the size at which the widths it has drawn would just meet the target
#' (the width of a Wald interval shrinks as \eqn{1 / \sqrt{N - 1}}{1 /
#' sqrt(N - 1)}), and from there up one case at a time to the first size
#' that meets it.
#'
#' Each sample is treated as the estimation functions treat it, so the
#' search counts the intervals a user would be given. A sample for which
#' \code{\link{reliability_omega}} or the model implied
#' \code{\link{reliability_alpha}} would return no interval (a
#' single-factor fit that does not converge, three items with no
#' single-factor solution, or a solution that has run away, whose error
#' variances sum below zero so that the model implied coefficient is above
#' one; McDonald, 2004) gave no interval, and so no narrow one: it is
#' counted as wider than \code{width}. So is a sample whose likelihood has
#' no maximum, whatever the sum, one error variance running off below
#' zero while the likelihood is flat or still rising, which the estimation
#' functions refuse with an error of class \code{dmar_no_maximum} (see
#' \code{\link{reliability_omega}}). The assurance is the probability
#' that the study yields an interval no wider than \code{width} (Terry &
#' Kelley, 2012, p. 373), not that probability among the studies that
#' yield an interval. When some samples of the final pass returned no
#' interval, the result says how many (the \code{no_interval} attribute,
#' printed under the table), since at the planned \emph{N} that share of
#' studies would end without a coefficient. Each confidence limit is clamped
#' into [0, 1], as the estimation functions clamp it, so a negative alpha
#' whose interval lies wholly below zero has the interval [0, 0]. A size
#' at which no sample returns an interval stops the search with the first
#' error.
#'
#' The cost of the assurance search depends on the interval. The normal
#' theory intervals are closed forms, so a search at the default
#' \code{G_initial = 500} and \code{G = 5000} finishes in
#' seconds. With the factor analytic interval a one factor model is fit
#' at every Monte Carlo iteration and at every candidate sample size the
#' search visits, so a congeneric plan with an assurance runs for tens of
#' seconds even at \code{G_initial = 50} and \code{G = 200}
#' and for many minutes at the defaults. The examples on this page
#' therefore stop at the closed form for the congeneric plan; the
#' assurance version is the same call with \code{assurance} supplied
#' (for example, \code{assurance = .80}). A \code{seed} makes the search
#' reproducible and leaves the caller's random number stream as it was;
#' without one the search draws from, and advances, the current stream.
#'
#' @section Planning Under Missing Data:
#' The fraction of missing information (FMI) of a parameter is
#' \deqn{\gamma = 1 - \frac{SE_C^2}{SE_O^2},}{gamma = 1 - SE_C^2 / SE_O^2,}
#' where \eqn{SE_O} is the standard error of the estimate from the data as
#' observed and \eqn{SE_C} the standard error the same sample would have
#' given had it been complete (Savalei & Rhemtulla, 2012, Eq. 4). It is
#' the proportion of the information about the parameter that the missing
#' values take away, one minus the asymptotic efficiency of the incomplete
#' data estimator relative to the complete data one. At a fixed \emph{N},
#' missing values therefore multiply the sampling variance of the estimate
#' by \eqn{1 / (1 - \gamma)} and the width of a Wald interval by
#' \eqn{1 / \sqrt{1 - \gamma}}{1 / sqrt(1 - gamma)}, the width inflation
#' factor of Savalei and Rhemtulla (2012); an FMI of .56 makes the interval
#' 1.5 times as wide.
#'
#' Every closed form on this page has a sampling variance of the form
#' \eqn{V / (N - 1)}, so its information about the coefficient is
#' \eqn{(N - 1) / V}. With missing data that information becomes
#' \eqn{(N - 1)(1 - \gamma) / V}, and the width target is met when
#' \eqn{(N - 1)(1 - \gamma)} equals \eqn{N_c - 1}, where \eqn{N_c} is the
#' complete data solution of the closed form before it is rounded up to a
#' whole number:
#' \deqn{N = \frac{N_c - 1}{1 - \gamma} + 1.}{N = (N_c - 1) / (1 - gamma) + 1.}
#' The factor scales \eqn{N - 1} rather than \emph{N} because \eqn{N - 1}
#' is what the planner's information is proportional to; the target sample
#' size of Savalei and Rhemtulla (2012), \eqn{N / (1 - \gamma)}, is the
#' same idea for an information proportional to \emph{N}, and the two
#' differ by \eqn{\gamma / (1 - \gamma)} cases (one case at an FMI of
#' .5). An FMI of .5 doubles \eqn{N - 1}, not \emph{N}. The plan is the
#' smallest whole \emph{N} at or above that value. The complete data sample
#' size, the whole number \eqn{N_c} rounds up to, travels with the result
#' as the \code{necessary_N_complete} attribute and is printed under the
#' table. Because the plan maps the unrounded \eqn{N_c}, the formula
#' applied to the printed whole number can give a larger \emph{N} than the
#' plan, by up to about \eqn{1 / (1 - \gamma)}{1 / (1 - gamma)} cases (the
#' population of the first example below needs 38 cases with complete data
#' and 42 at an FMI of .1, while the formula applied to 38 gives 43).
#'
#' With an \code{assurance}, the Monte Carlo search runs exactly as it does
#' for complete data, starting from the complete data closed form, and the
#' size \eqn{N_e} it finds is mapped back with the same formula. The search
#' therefore takes the distribution of the interval width at \emph{N} under
#' missing data to be the complete data width distribution at the effective
#' size \eqn{(N - 1)(1 - \gamma) + 1}, which keeps the larger relative
#' variability of the widths of a smaller sample (scaling the complete data
#' widths at \emph{N} by the width inflation factor would not). This is a
#' first-order approximation, and
#' \code{\link{ss_aipe_reliability_sensitivity}} checks it by simulating
#' incomplete data and analyzing them by full information maximum
#' likelihood (its \code{fmi} and \code{missing_prop} arguments). The
#' \code{assurance_attained} row is the assurance found at \eqn{N_e}, the
#' \code{necessary_N_complete} attribute is \eqn{N_e} itself (a whole
#' number, so the formula applied to it reproduces the plan), and a
#' \code{start_N} is read as a size under missing data and converted to
#' its effective size before the search.
#'
#' \strong{Where the FMI comes from.} The FMI belongs to the coefficient,
#' not to the items, and it is not the proportion of missing values: it
#' depends on the missing data mechanism, the model, and how strongly the
#' items are related. Information lost on one item is partly recovered
#' from the items correlated with it, so the FMI can be lower than the rate
#' of missing data; under data missing completely at random it is bounded
#' by the largest rates of univariate and pairwise missingness, and under
#' data missing at random it can be well above the rate of missing data
#' (Savalei & Rhemtulla, 2012). At the design stage it comes from a pilot
#' study or a comparable data set analyzed the way the planned study will
#' be: fit the single-factor model by full information maximum likelihood
#' (for example \code{\link{fa_confirmatory}(..., missing = "fiml", output = "fit")},
#' with \code{equal_loading = TRUE} for the tau equivalent coefficient that
#' the Parallel and True Score plans target) and take the FMI of its
#' \code{omega_} defined parameter from
#' \code{\link{fraction_missing_information}}. That FMI refers to the cases
#' the fit used: a case with no observed item carries no information, and
#' the fit sets it aside (the \code{N_excluded} attribute counts such
#' cases). When the planned \emph{N} counts every case recruited, those
#' who will answer no item included, fold them in as
#' \eqn{1 - (N_{used} / N)(1 - \gamma)}{1 - (N_used / N)(1 - gamma)}, with
#' \eqn{N_{used}}{N_used} the cases the fit used and \emph{N} those
#' recruited; \code{\link{ss_aipe_reliability_sensitivity}} reports the
#' FMI its simulated data realized in the same way. A sample FMI is
#' imprecise. Chen
#' and Savalei (2021) found it close to unbiased even at \eqn{N = 50} but
#' variable from sample to sample: precise estimates took several hundred
#' cases at overall missing rates of 10 to 20 percent and more than 1,000
#' at 30 percent. Planning with the upper limit of an interval for the
#' pilot's FMI, rather than with its point estimate, is the conservative
#' choice. Von Hippel (2020) recommends the same use of a pilot's upper
#' limit when deciding how many imputations to run; carrying it over to
#' sample size planning is our extension.
#'
#' Two things \code{fmi} leaves alone. The \code{data} argument still
#' estimates the planning covariance matrix from the complete cases of the
#' pilot data
#' (\code{na.omit()}); \code{fmi} is what carries the missing data of the
#' planned study. And the planned interval is the same Wald interval with a
#' normal critical value; the FMI changes how much information the sample
#' carries, not the reference distribution.
#'
#' @param model The measurement model assumed for the population.
#'   Accepts (case-sensitive aliases shown in parentheses):
#'   \code{"Parallel"} (\code{"parallel"}, \code{"SB"},
#'   \code{"Spearman Brown"}, \code{"Spearman-Brown"}, \code{"sb"}) for
#'   the strictly parallel items model; \code{"True Score"}
#'   (\code{"True Score Equivalent"}, \code{"True-Score Equivalent"},
#'   \code{"Equivalent"}, \code{"Tau Equivalent"}, \code{"tau-equivalent"},
#'   \code{"Tau-Equivalent"}, \code{"True-Score"}, \code{"true-score"},
#'   \code{"true score"}, \code{"Cronbach"}, \code{"cronbach"},
#'   \code{"Chronbach"}, \code{"alpha"}) for the tau equivalent model
#'   (in which case the function plans for coefficient alpha); or
#'   \code{"Congeneric"} (\code{"congeneric"}, \code{"omega"},
#'   \code{"Omega"}) for the congeneric model (in which case the function
#'   plans for coefficient omega).
#' @param type The method used to construct the confidence interval on the
#'   reliability coefficient: either \code{"Factor Analytic"} (McDonald,
#'   1999), available for all three measurement models, or
#'   \code{"Normal Theory"} (van Zyl, Neudecker, & Nel, 2000), available
#'   for the parallel and tau equivalent models only (with the congeneric
#'   model it is an error).
#' @param width The desired full width of the two-sided confidence interval.
#' @param true_cov_matrix A symmetric population covariance matrix among
#'   the items (a planning value, not a sample covariance matrix), used
#'   as it is to imply the population reliability and its sampling
#'   distribution; a matrix with a negative eigenvalue is refused. A
#'   correlation matrix is a covariance matrix of standardized items and
#'   plans for them.
#' @param conf_level Confidence level (i.e., \eqn{1 - \alpha}, where
#'   \eqn{\alpha} is the Type I error rate). Default \code{0.95}.
#' @param assurance Optional probability with which the realized interval
#'   is to be no wider than \code{width} (a sample that returns no interval
#'   counts as one that is wider). When \code{NULL} (the default),
#'   the planner targets the \emph{expected} width; when supplied
#'   (e.g., 0.80, 0.85, 0.95), the function follows the closed form with a
#'   Monte Carlo search.
#' @param data A data set from which the population covariance matrix
#'   should be inferred, from its complete cases; the number of items is
#'   the number of its columns.
#' @param n_items Number of items (\emph{J} in the formulas of the Details).
#'   Needed when the population is given by \code{true_item_covariance} or
#'   \code{true_factor_loading} with \code{true_residual_var}; with
#'   \code{true_cov_matrix} or \code{data} it is the number of columns and
#'   need not be given.
#' @param true_item_covariance The presumed common covariance \emph{c} of
#'   every pair of items, which with \code{true_residual_var} gives the
#'   population: the common loading is the square root of \emph{c}, so two
#'   items covary \emph{c} and an item's variance is \emph{c} plus its
#'   \code{true_residual_var}. It is the inter-item correlation of Terry
#'   and Kelley (2012) when the items are standardized, that is, when
#'   \code{true_residual_var} is \eqn{1 - c}; otherwise items \emph{j} and
#'   \emph{k} correlate
#'   \eqn{c / \sqrt{(c + \psi^2_j)(c + \psi^2_k)}}{c / sqrt((c + psi_j^2)(c
#'   + psi_k^2))}, with \eqn{\psi^2_j} the residual variance of item
#'   \emph{j}. With \code{true_item_covariance = .3} and
#'   \code{true_residual_var = .2}, for example, the items correlate .6, not
#'   .3. One value.
#' @param true_factor_loading Population factor loading or loadings: one
#'   number for the Parallel and True Score models, whose items share one
#'   loading, and a vector of \code{n_items} loadings for the Congeneric
#'   model. A vector for the Parallel or True Score model is an error that
#'   names the Congeneric model, since unequal loadings are what
#'   distinguish it.
#' @param true_residual_var Population residual (unique, or error) variance
#'   or variances: one number for the Parallel model, whose items share
#'   one, and a vector of \code{n_items} variances for the True Score and
#'   Congeneric models, each
#'   greater than zero; with the normal theory interval of the True Score
#'   model, any values that leave the covariance matrix they imply a
#'   covariance matrix (see Details).
#' @param G_initial Number of Monte Carlo replications in each pass of the
#'   initial assurance search (default 500): a single whole number of at
#'   least 1.
#' @param G Number of Monte Carlo replications in each pass of the final
#'   assurance verification (default 5000), the count the reported
#'   \code{assurance_attained} rests on: a single whole number of at
#'   least 1. Anything else given for either count is refused with an error
#'   of class \code{dmar_count_invalid}. Neither is used without
#'   \code{assurance}.
#' @param start_N Optional starting total sample size for the iterative
#'   assurance search.
#' @param verbose If \code{TRUE}, prints the current sample size and
#'   empirical assurance at each step of the Monte Carlo search.
#' @param fmi The fraction of missing information anticipated for the
#'   reliability coefficient when it is estimated by full information
#'   maximum likelihood from incomplete data: a single number in
#'   \eqn{[0, 1)}, with the default \code{0} for complete data, which
#'   returns exactly the complete data plan. It is not the proportion of
#'   missing values. Writing \eqn{\gamma} for \code{fmi}, the plan
#'   divides the complete data \eqn{N - 1} by \eqn{1 - \gamma}; see the
#'   section on planning under missing data for the derivation and for
#'   where a planning value comes from.
#' @param seed Optional seed for the Monte Carlo assurance search. The
#'   default \code{NULL} draws from the current state of the random number
#'   generator; a number makes the search reproducible, and the caller's
#'   generator state is restored when the function returns. The closed
#'   form (no \code{assurance}) draws nothing.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}.
#' Without \code{assurance} the data frame has a single row,
#' \code{"necessary_N"}, giving the necessary \emph{N}. With \code{assurance}
#' supplied, the data frame has five rows: \code{"necessary_N"} (necessary
#' \emph{N}), \code{"width"} (echo of the target width),
#' \code{"assurance"} (echo of the requested probability),
#' \code{"assurance_attained"} (the assurance achieved at the returned
#' \emph{N} in the Monte Carlo verification), and \code{"G"}
#' (the number of Monte Carlo replications of the final pass). When \code{fmi} is above 0,
#' \code{necessary_N} is the size under missing data, the rows are
#' otherwise unchanged, and two attributes record the plan:
#' \code{fmi}, the value planned for, and \code{necessary_N_complete},
#' the size complete data would need (for the closed form, its solution
#' rounded up, while the plan maps the unrounded solution; for an assurance
#' plan, the size at which the assurance was simulated, so
#' \code{assurance_attained} refers to it). With \code{fmi = 0} neither
#' attribute is added. An assurance plan carries the number of samples of
#' its final Monte Carlo pass that returned no interval as the
#' \code{no_interval} attribute (zero when every sample returned one),
#' printed under the table when above zero; they are counted as wide in
#' \code{assurance_attained}. Both forms carry the class
#' \code{dmar_ss_aipe} of the AIPE planners: \code{generics::tidy()} gives
#' the planned \emph{N}, beside the target width for an assurance plan
#' (the closed form's table holds only \code{necessary_N}, so its width is
#' \code{NA}), and \code{generics::glance()} adds the
#' remaining rows as columns, and the \code{fmi},
#' \code{necessary_N_complete}, and \code{no_interval} attributes when the
#' plan has them (see \code{\link{dmar_tidiers}}).
#'
#' @references
#' Chen, L., & Savalei, V. (2021). Three sample estimates of fraction of
#'   missing information from full information maximum likelihood.
#'   \emph{Frontiers in Psychology, 12}, 667802.
#'   \doi{10.3389/fpsyg.2021.667802}
#'
#' Cronbach, L. J. (1951). Coefficient alpha and the internal structure
#'   of tests. \emph{Psychometrika, 16}(3), 297--334.
#'
#' Guttman, L. (1945). A basis for analyzing test-retest reliability.
#'   \emph{Psychometrika, 10}(4), 255--282. \doi{10.1007/BF02288892}
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
#' McDonald, R. P. (1999). \emph{Test theory: A unified treatment}. Mahwah, NJ:
#'   Lawrence Erlbaum Associates.
#'
#' McDonald, R. P. (2004). Respecifying improper structures.
#'   \emph{Structural Equation Modeling, 11}(2), 194--209.
#'
#' Raykov, T. (2002). Analytic estimation of standard error and confidence
#'   interval for scale reliability.
#'   \emph{Multivariate Behavioral Research, 37}, 89--103.
#'   \doi{10.1207/S15327906MBR3701_04}
#'
#' Savalei, V., & Rhemtulla, M. (2012). On obtaining estimates of the
#'   fraction of missing information from full information maximum
#'   likelihood. \emph{Structural Equation Modeling, 19}(3), 477--494.
#'   \doi{10.1080/10705511.2012.687669}
#'
#' Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
#'   reliability coefficients: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{British Journal of Mathematical
#'   and Statistical Psychology, 65}, 371--401.
#'   \doi{10.1111/j.2044-8317.2011.02030.x}
#'
#' van Zyl, J. M., Neudecker, H., & Nel, D. G. (2000). On the distribution of
#'   the maximum likelihood estimator of Cronbach's alpha.
#'   \emph{Psychometrika, 65}(3), 271--280. \doi{10.1007/BF02296146}
#'
#' von Hippel, P. T. (2020). How many imputations do you need? A two-stage
#'   calculation using a quadratic rule. \emph{Sociological Methods &
#'   Research, 49}(3), 699--718. \doi{10.1177/0049124117747303}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note Not all of the items can be entered into the function to represent the population values.
#' For example, either \code{data} can be used, or \code{true_cov_matrix},
#' or \code{n_items}, \code{true_item_covariance}, and
#' \code{true_residual_var}, or \code{n_items}, \code{true_factor_loading},
#' and \code{true_residual_var}. With a large number of replications
#' (\code{G}) this function may take considerable time.
#'
#' @section Warning: In some conditions the factor analytic model fit by
#' \code{\link{cfa_1}} (via \pkg{lavaan}) may fail to converge, and you may
#' see a non-convergence message from \pkg{lavaan}. The Monte Carlo assurance
#' search counts a non-converged iteration as an interval wider than
#' \code{width} and continues, and the number of such iterations in the
#' final pass is printed under the result; a sample size at which no
#' iteration converges stops the search with the first error. Frequent
#' non-convergence usually means the model is poorly determined by the
#' data, for example because of a small sample size or a poorly behaved
#' covariance matrix, and it raises the planned \emph{N}, since a study
#' that yields no interval has not yielded a narrow one.
#'
#' @seealso \code{\link{reliability}}, \code{\link{cfa_1}},
#'   \code{\link{ss_aipe_reliability_sensitivity}} to check a plan, under
#'   complete or missing data, by simulation.
#'
#' @examples
#' # Expected confidence interval width (closed form, no Monte Carlo search)
#' # for six parallel items that correlate .6: standardized items, so
#' # true_residual_var is 1 - true_item_covariance.
#' ss_aipe_reliability(model = "Parallel", type = "Normal Theory", width = .1,
#'   n_items = 6, true_item_covariance = .6, true_residual_var = .4,
#'   conf_level = .95, assurance = NULL)
#'
#' # The same plan when the coefficient will be estimated by full information
#' # maximum likelihood and is expected to lose 30 percent of its information
#' # to missing values (an FMI of .30 from a pilot study). The complete data
#' # N minus 1 is divided by 1 minus .30, and the complete data answer prints
#' # under the table.
#' ss_aipe_reliability(model = "Parallel", type = "Normal Theory", width = .1,
#'   n_items = 6, true_item_covariance = .6, true_residual_var = .4,
#'   conf_level = .95, fmi = .30)
#'
#' # The assurance cases run a Monte Carlo search. G_initial = 20 and
#' # G = 100 keep the examples quick; a reported plan deserves the
#' # defaults of 500 and 5000. The seed makes each search reproducible and
#' # leaves the random number stream of the session where it was.
#'
#' # Same population, now targeting an assurance.
#' ss_aipe_reliability(model = "Parallel", type = "Normal Theory", width = .1,
#'   n_items = 6, true_item_covariance = .6, true_residual_var = .4,
#'   conf_level = .95, assurance = .85, G_initial = 20, G = 100, seed = 113)
#'
#' # The assurance plan under the same FMI. The search runs at the complete
#' # data size, the necessary N of the call above (kept as the attribute
#' # necessary_N_complete), and maps it back.
#' ss_aipe_reliability(model = "Parallel", type = "Normal Theory", width = .1,
#'   n_items = 6, true_item_covariance = .6, true_residual_var = .4,
#'   conf_level = .95, assurance = .85, G_initial = 20, G = 100, seed = 113,
#'   fmi = .30)
#'
#' # The true score (tau equivalent) model takes true_residual_var as a
#' # vector of length n_items (number of items) while true_item_covariance
#' # stays a single value. Error
#' # variances around 1 - .5 give item variances from .9 to 1.1 and
#' # correlations near .5, and the plan is for alpha of these items as they
#' # are, unequal variances and all.
#' ss_aipe_reliability(model = "True Score", type = "Normal Theory",
#'   width = .1, n_items = 5, true_item_covariance = .5,
#'   true_residual_var = c(.4, .5, .6, .4, .6),
#'   conf_level = .95, assurance = .85, G_initial = 20, G = 100,
#'   seed = 113)
#'
#' # Congeneric model, planned from the item loadings and error variances rather
#' # than from a single covariance. With assurance = NULL the necessary N comes
#' # from the delta method variance of the model implied omega at the
#' # population, the variance behind reliability_omega(denominator =
#' # "model_implied", ci_method = "ml"). Adding an assurance to this plan fits
#' # a one factor model at every Monte Carlo iteration and is the slow case,
#' # so the page stops at the closed form (see Details).
#' ss_aipe_reliability(model = "Congeneric", type = "Factor Analytic", width = .15,
#'   n_items = 4, true_factor_loading = c(.8, .7, .7, .8),
#'   true_residual_var = c(.4, .5, .5, .4),
#'   conf_level = .95, assurance = NULL)
#'
#' # Planning from a presumed population correlation matrix among the items,
#' # which plans for standardized items.
#' pop_mat <- rbind(
#'   c(1.0000000, 0.3813850, 0.4216370, 0.3651484, 0.4472136),
#'   c(0.3813850, 1.0000000, 0.4020151, 0.3481553, 0.4264014),
#'   c(0.4216370, 0.4020151, 1.0000000, 0.3849002, 0.4714045),
#'   c(0.3651484, 0.3481553, 0.3849002, 1.0000000, 0.4082483),
#'   c(0.4472136, 0.4264014, 0.4714045, 0.4082483, 1.0000000))
#' ss_aipe_reliability(model = "True Score", type = "Normal Theory", width = .15,
#'   true_cov_matrix = pop_mat, conf_level = .95, assurance = .85,
#'   G_initial = 20, G = 100, seed = 113)
#'
#' @keywords design
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export

ss_aipe_reliability <- function(model = NULL, type = NULL, width = NULL,
                                true_cov_matrix = NULL,
                                conf_level = 0.95, assurance = NULL, data = NULL, n_items = NULL,
                                true_item_covariance = NULL,
                                true_factor_loading = NULL,
                                true_residual_var = NULL,
                                G_initial = 500, G = 5000, start_N = NULL, verbose = FALSE,
                                fmi = 0, seed = NULL) {
  .dmar_check_conf_level(conf_level)
  if (!is.null(width)) .dmar_check_width(width)
  .dmar_check_assurance(assurance)
  .dmar_check_numbers("ss_aipe_reliability", c("n_items", "true_item_covariance"))
  if (!requireNamespace("MASS", quietly = TRUE)) .dmar_refuse(
    "The package 'MASS' is needed; please install the package and try again.",
    "package_required", package = "MASS")
  .aipe_check_fmi(fmi)
  # The two counts of the assurance search are checked whether or not an
  # assurance is asked for, as every count of replications in the package is.
  G_initial <- .dmar_check_count(G_initial, "G_initial")
  G <- .dmar_check_count(G, "G")
  # The closed forms draw nothing; the assurance search below is where the
  # random draws begin, and a supplied seed is restored on exit so the
  # caller's stream is left as it was.
  .dmar_local_seed(seed)

  model_type1 <- c(
    "Parallel", "SB", "Spearman Brown", "Spearman-Brown",
    "sb", "parallel"
  )
  model_type2 <- c(
    "True Score", "True Score Equivalent", "True-Score Equivalent",
    "Equivalent", "Tau Equivalent", "Chronbach", "Cronbach",
    "cronbach", "alpha", "true score", "tau-equivalent",
    "Tau-Equivalent", "True-Score", "true-score"
  )
  model_type3 <- c("Congeneric", "congeneric", "omega", "Omega")
  if (sum(model == model_type1, model == model_type2, model ==
    model_type3) != 1) {
    .dmar_refuse_argument(
      "Assign one and only one of the three types of models to 'model'.",
      arg = "model")
  }
  if (sum(model == model_type1) == 1) {
    Model_to_Use <- "Parallel"
  }
  if (sum(model == model_type2) == 1) {
    Model_to_Use <- "True Score"
  }
  if (sum(model == model_type3) == 1) {
    Model_to_Use <- "Congeneric"
  }
  type_1 <- c(
    "Normal Theory", "Normal theory", "normal theory",
    "nt", "NT"
  )
  type_2 <- c(
    "Factor Analytic", "factor analytic", "Factor analytic",
    "factor Analytic"
  )
  if (sum(type == type_1, type == type_2) != 1) {
    .dmar_refuse_argument(
      "Assign either Factor Analytic or Normal Theory to 'type'.", arg = "type")
  }
  if (sum(type == type_1) == 1) {
    Type_to_Use <- "Normal Theory"
  }
  if (sum(type == type_2) == 1) {
    Type_to_Use <- "Factor Analytic"
  }
  # The congeneric model has no normal theory interval. The refusal comes
  # before anything is computed: the closed form below does not read type,
  # and in the assurance search every iteration's interval stopped inside
  # try(), so the search failed with an unrelated error.
  if (Model_to_Use == "Congeneric" && Type_to_Use == "Normal Theory") {
    .dmar_refuse_argument(c(
      "The congeneric model has no normal theory interval: type = ",
      "\"Normal Theory\" is for the Parallel and True Score models, whose ",
      "coefficient is alpha. Use type = \"Factor Analytic\" to plan for ",
      "coefficient omega."), arg = NULL)
  }
  if (sum(model == model_type1) == 1) {
    if (!is.null(true_item_covariance)) {
      if (!is.null(true_factor_loading)) {
        .dmar_refuse_argument(
          "Give 'true_item_covariance' or 'true_factor_loading', not both.",
          arg = c("true_item_covariance", "true_factor_loading"))
      }
      if (is.null(true_residual_var)) {
        .dmar_refuse_argument(
          "'n_items', 'true_residual_var', and one of 'true_item_covariance' and 'true_factor_loading' must all be given.",
          arg = c("n_items", "true_residual_var", "true_item_covariance", "true_factor_loading"))
      }
      # The length is checked before the value: a vector true_residual_var would
      # otherwise stop in the comparison below with R's own "the condition
      # has length > 1" rather than with this message.
      if (length(true_residual_var) > 1) {
        .dmar_refuse_argument(
          "'true_residual_var' must be one number for the Parallel model, whose items have one error variance. For error variances that differ, use model = \"True Score\" or model = \"Congeneric\".",
          arg = "true_residual_var")
      }
      if (true_residual_var <= 0) {
        .dmar_refuse_argument(
          "'true_residual_var' must be greater than zero.", arg = "true_residual_var")
      }
      if (is.null(n_items)) {
        .dmar_refuse_argument(
          "'n_items', 'true_residual_var', and one of 'true_item_covariance' and 'true_factor_loading' must all be given.",
          arg = c("n_items", "true_residual_var", "true_item_covariance", "true_factor_loading"))
      }
      if (n_items <= 0) {
        .dmar_refuse_argument(
          "'n_items' must be greater than zero.", arg = "n_items")
      }
      if (length(true_item_covariance) >= 2) {
        .dmar_refuse_argument(
          "'true_item_covariance' must be a single number.",
          arg = "true_item_covariance")
      }
      lambda_1 <- sqrt(true_item_covariance)
      v_lambda <- matrix(data = lambda_1, nrow = 1, ncol = n_items)
      v_Psi <- matrix(data = true_residual_var, nrow = 1, ncol = n_items)
      temp_mat <- covmat_from_cfa(
        factor_loading = v_lambda, residual_var = v_Psi,
        tol_det = 1e-05
      )
      cov_mat <- Sigma <- temp_mat$model_implied_cov
    }
    if (!is.null(true_factor_loading)) {
      if (!is.null(true_item_covariance)) {
        .dmar_refuse_argument(
          "Give 'true_item_covariance' or 'true_factor_loading', not both.",
          arg = c("true_item_covariance", "true_factor_loading"))
      }
      # Parallel items share one loading. A vector was recycled across the
      # items by matrix() below, which planned a congeneric loading pattern
      # under the Parallel label without saying so.
      if (length(true_factor_loading) != 1) {
        .dmar_refuse_argument(
          "'true_factor_loading' must be one number for the Parallel model. If you want to enter multiple 'true_factor_loading' values, you must use the Congeneric model.",
          arg = "true_factor_loading")
      }
      if (is.null(true_residual_var)) {
        .dmar_refuse_argument(
          "'n_items', 'true_residual_var', and one of 'true_item_covariance' and 'true_factor_loading' must all be given.",
          arg = c("n_items", "true_residual_var", "true_item_covariance", "true_factor_loading"))
      }
      if (length(true_residual_var) > 1) {
        .dmar_refuse_argument(
          "'true_residual_var' must be one number for the Parallel model, whose items have one error variance. For error variances that differ, use model = \"True Score\" or model = \"Congeneric\".",
          arg = "true_residual_var")
      }
      if (true_residual_var <= 0) {
        .dmar_refuse_argument(
          "'true_residual_var' must be greater than zero.", arg = "true_residual_var")
      }
      if (is.null(n_items)) {
        .dmar_refuse_argument(
          "'n_items', 'true_residual_var', and one of 'true_item_covariance' and 'true_factor_loading' must all be given.",
          arg = c("n_items", "true_residual_var", "true_item_covariance", "true_factor_loading"))
      }
      if (n_items <= 0) {
        .dmar_refuse_argument(
          "'n_items' must be greater than zero.", arg = "n_items")
      }
      v_lambda <- matrix(data = true_factor_loading, nrow = 1, ncol = n_items)
      v_Psi <- matrix(data = true_residual_var, nrow = 1, ncol = n_items)
      temp_mat <- covmat_from_cfa(
        factor_loading = v_lambda, residual_var = v_Psi,
        tol_det = 1e-05
      )
      cov_mat <- Sigma <- temp_mat$model_implied_cov
    }
    if (!is.null(true_cov_matrix)) {
      if (!isSymmetric(true_cov_matrix, tol = 1e-05)) {
        .dmar_refuse_argument(
          "Input a symmetric covariance or correlation matrix, 'true_cov_matrix'",
          arg = "true_cov_matrix")
      }
      .ss_aipe_check_S(true_cov_matrix)
      cov_mat <- true_cov_matrix
      n_items <- dim(true_cov_matrix)[1]
    }
    if (!is.null(data)) {
      data <- na.omit(data)
      cov_mat <- var(data, y = NULL, na.rm = TRUE)
    }
    if (!is.null(true_cov_matrix) || !is.null(data)) {
      .ss_aipe_check_parallel(cov_mat)
    }
    sigma_jj <- sum(diag(cov_mat))
    sigma2_Y <- sum(cov_mat)
    p <- ncol(cov_mat)
    alpha <- (p / (p - 1)) * (1 - sigma_jj / sigma2_Y)
    # The factor analytic interval is that of the coefficient implied by
    # the single-factor model with equal loadings, which is zero for items
    # that share no common factor, and the closed form below is a function
    # of alpha; see .ss_aipe_check_common_factor().
    if (Type_to_Use == "Factor Analytic") {
      .ss_aipe_check_common_factor(cov_mat, parallel = TRUE,
                                   fitted = !is.null(true_cov_matrix) || !is.null(data))
    }
    k <- 1 - alpha
    Crit_Value <- qnorm(1 - (1 - conf_level) / 2)
    Nec_N <- as.numeric(((Crit_Value^2) * 8 * (k^2) * p) / ((width^2) *
      (p - 1)) + 1)
  }
  if (sum(model == model_type2) == 1) {
    # The equal loading and the error variances of the population, when it
    # is given by them; a population given by true_cov_matrix or data is fit when the
    # factor analytic closed form below needs them.
    pop_lambda <- pop_psi <- NULL
    if (!is.null(true_item_covariance)) {
      if (!is.null(true_factor_loading)) {
        .dmar_refuse_argument(
          "Give 'true_item_covariance' or 'true_factor_loading', not both.",
          arg = c("true_item_covariance", "true_factor_loading"))
      }
      if (is.null(true_residual_var)) {
        .dmar_refuse_argument(
          "'true_residual_var' must be given: the error variance of each item.",
          arg = "true_residual_var")
      }
      if (is.null(n_items)) {
        .dmar_refuse_argument(
          "'n_items', the number of items, must be given.", arg = "n_items")
      }
      if ((n_items) != length(true_residual_var)) {
        .dmar_refuse_argument(
          "'true_residual_var' must have 'n_items' values, one for each item.",
          arg = c("true_residual_var", "n_items"))
      }
      if (n_items <= 0) {
        .dmar_refuse_argument(
          "'n_items' must be greater than zero.", arg = "n_items")
      }
      if (length(true_item_covariance) >= 2) {
        .dmar_refuse_argument(
          "You can only enter one 'true_item_covariance' value.",
          arg = "true_item_covariance")
      }
      lambda_1 <- sqrt(true_item_covariance)
      lambda_vector <- rep(lambda_1, times = n_items)
      .ss_aipe_check_psi(true_residual_var, lambda_vector,
                         normal_theory = Type_to_Use == "Normal Theory")
      Population_Cov <- covmat_from_cfa(
        factor_loading = lambda_vector,
        residual_var = true_residual_var
      )$model_implied_cov
      cov_mat <- Population_Cov
      pop_lambda <- lambda_vector
      pop_psi <- true_residual_var
    }
    if (!is.null(true_factor_loading)) {
      if (!is.null(true_item_covariance)) {
        .dmar_refuse_argument(
          "Give 'true_item_covariance' or 'true_factor_loading', not both.",
          arg = c("true_item_covariance", "true_factor_loading"))
      }
      if (is.null(true_residual_var)) {
        .dmar_refuse_argument(
          "'n_items', 'true_residual_var', and one of 'true_item_covariance' and 'true_factor_loading' must all be given.",
          arg = c("n_items", "true_residual_var", "true_item_covariance", "true_factor_loading"))
      }
      # n_items is checked for presence before it is compared with the length of
      # true_residual_var: with n_items NULL the comparison is empty and stops with R's
      # own "argument is of length zero" rather than with this message.
      if (is.null(n_items)) {
        .dmar_refuse_argument(
          "You need to enter a quantity for n_items if you enter true_factor_loading and true_residual_var",
          arg = NULL)
      }
      if (length(true_residual_var) != (n_items)) {
        .dmar_refuse_argument(
          "'true_residual_var' must have 'n_items' values, one for each item.",
          arg = c("true_residual_var", "n_items"))
      }
      if (n_items <= 0) {
        .dmar_refuse_argument("'n_items' must be greater than zero.", arg = "n_items")
      }
      if (length(true_factor_loading) != 1) {
        .dmar_refuse_argument(
          "'true_factor_loading' must be one number for the True Score model. If you want to enter multiple 'true_factor_loading' values, you must use the Congeneric model.",
          arg = "true_factor_loading")
      }
      lambda_vector <- rep(true_factor_loading, times = n_items)
      .ss_aipe_check_psi(true_residual_var, lambda_vector,
                         normal_theory = Type_to_Use == "Normal Theory")
      Population_Cov <- covmat_from_cfa(
        factor_loading = lambda_vector,
        residual_var = true_residual_var
      )$model_implied_cov
      cov_mat <- Population_Cov
      pop_lambda <- lambda_vector
      pop_psi <- true_residual_var
    }
    if (!is.null(true_cov_matrix)) {
      if (!isSymmetric(true_cov_matrix, tol = 1e-05)) {
        .dmar_refuse_argument(
          "Input a symmetric covariance or correlation matrix 'true_cov_matrix'",
          arg = "true_cov_matrix")
      }
      .ss_aipe_check_S(true_cov_matrix)
      cov_mat <- true_cov_matrix
      n_items <- dim(true_cov_matrix)[1]
      pop_lambda <- pop_psi <- NULL
    }
    if (!is.null(data)) {
      data <- na.omit(data)
      cov_mat <- var(data, y = NULL, na.rm = TRUE)
      pop_lambda <- pop_psi <- NULL
    }
    Crit_Value <- qnorm(1 - (1 - conf_level) / 2)
    if (Type_to_Use == "Normal Theory") {
      # Terry and Kelley (2012, Eq. 19): the van Zyl, Neudecker, and Nel
      # (2000) variance of alpha, evaluated on the population covariance
      # matrix itself (Eq. 8), so the plan is for the alpha that
      # reliability_alpha() estimates from the raw item scores. Before
      # 1.1.0 the matrix was first turned into a correlation matrix, which
      # planned for the alpha of standardized items.
      V <- .alpha_vanzyl_avar(cov_mat)
    } else {
      # The factor analytic interval of the tau equivalent model is the
      # delta method Wald interval of the coefficient implied by the
      # single-factor model with equal loadings, the interval of
      # reliability_alpha(estimator = "model_implied", ci_method = "ml")
      # and the one the assurance search draws, so its variance is the one
      # planned with (.ss_aipe_omega_avar(equal_loadings = TRUE)). The
      # model implied coefficient is the maximum likelihood estimator under
      # the model, so its variance is no larger than the normal theory
      # variance of alpha, and smaller whenever the error variances differ:
      # for loadings .7 and error variances .1, .2, 1, and 2 at a width of
      # .10 the plan is 242 cases, against 367 for the normal theory
      # interval on the covariance matrix. With equal error variances the
      # two coincide. Before 1.1.0 the closed form was the normal theory
      # one whatever the type, evaluated on the correlation matrix, which
      # gave 206 cases there: too few for the factor analytic interval of
      # these items, which averages .108 wide at 206 cases. A population
      # given by true_cov_matrix or data is first fit by the model.
      pop <- .ss_aipe_check_common_factor(cov_mat, pop_lambda, pop_psi,
                                          fitted = is.null(pop_lambda))
      V <- .ss_aipe_omega_avar(cov_mat, pop$lambda, pop$psi,
                               equal_loadings = TRUE)
    }
    Nec_N <- V / (width / (2 * Crit_Value))^2 + 1
  }
  if (sum(model == model_type3) == 1) {
    if (!is.null(true_item_covariance)) {
      if (!is.null(true_factor_loading)) {
        .dmar_refuse_argument(
          "Give 'true_item_covariance' or 'true_factor_loading', not both.",
          arg = c("true_item_covariance", "true_factor_loading"))
      }
      if (is.null(true_residual_var)) {
        .dmar_refuse_argument(
          "'n_items', 'true_residual_var', and one of 'true_item_covariance' and 'true_factor_loading' must all be given.",
          arg = c("n_items", "true_residual_var", "true_item_covariance", "true_factor_loading"))
      }
      if (is.null(n_items)) {
        .dmar_refuse_argument(
          "'n_items', 'true_residual_var', and one of 'true_item_covariance' and 'true_factor_loading' must all be given.",
          arg = c("n_items", "true_residual_var", "true_item_covariance", "true_factor_loading"))
      }
      if (n_items <= 0) {
        .dmar_refuse_argument("'n_items' must be greater than zero.", arg = "n_items")
      }
      if ((n_items) != length(true_residual_var)) {
        .dmar_refuse_argument(
          "'true_residual_var' must have 'n_items' values, one for each item.",
          arg = c("true_residual_var", "n_items"))
      }
      if (length(true_item_covariance) >= 2) {
        .dmar_refuse_argument(
          "If you have multiple values for 'true_item_covariance', please put them as 'true_factor_loading' values. The square root of 'true_item_covariance' is the common loading.",
          arg = c("true_item_covariance", "true_factor_loading"))
      }
      .ss_aipe_check_psi(true_residual_var)
      lambda_1 <- sqrt(true_item_covariance)
      v_lambda <- matrix(data = lambda_1, nrow = 1, ncol = n_items)
      v_Psi <- matrix(data = true_residual_var, nrow = 1, ncol = n_items)
      temp_mat <- covmat_from_cfa(
        factor_loading = v_lambda, residual_var = v_Psi,
        tol_det = 1e-05
      )
      cov_mat <- Sigma <- temp_mat$model_implied_cov
      pop_lambda <- as.numeric(v_lambda)
      pop_psi <- as.numeric(v_Psi)
    }
    if (!is.null(true_factor_loading)) {
      if (!is.null(true_item_covariance)) {
        .dmar_refuse_argument(
          "Give 'true_item_covariance' or 'true_factor_loading', not both.",
          arg = c("true_item_covariance", "true_factor_loading"))
      }
      if (is.null(true_residual_var)) {
        .dmar_refuse_argument(
          "'n_items', 'true_residual_var', and one of 'true_item_covariance' and 'true_factor_loading' must all be given.",
          arg = c("n_items", "true_residual_var", "true_item_covariance", "true_factor_loading"))
      }
      if (is.null(n_items)) {
        .dmar_refuse_argument(
          "'n_items', 'true_residual_var', and one of 'true_item_covariance' and 'true_factor_loading' must all be given.",
          arg = c("n_items", "true_residual_var", "true_item_covariance", "true_factor_loading"))
      }
      if (n_items <= 0) {
        .dmar_refuse_argument("'n_items' must be greater than zero.", arg = "n_items")
      }
      if ((n_items) != length(true_factor_loading)) {
        .dmar_refuse_argument(
          "'true_factor_loading' must have 'n_items' values, one for each item.",
          arg = c("true_factor_loading", "n_items"))
      }
      if ((n_items) != length(true_residual_var)) {
        .dmar_refuse_argument(
          "'true_residual_var' must have 'n_items' values, one for each item.",
          arg = c("true_residual_var", "n_items"))
      }
      .ss_aipe_check_psi(true_residual_var)
      v_lambda <- matrix(data = true_factor_loading, nrow = 1, ncol = n_items)
      v_Psi <- matrix(data = true_residual_var, nrow = 1, ncol = n_items)
      temp_mat <- covmat_from_cfa(
        factor_loading = v_lambda, residual_var = v_Psi,
        tol_det = 1e-05
      )
      cov_mat <- Sigma <- temp_mat$model_implied_cov
      pop_lambda <- as.numeric(v_lambda)
      pop_psi <- as.numeric(v_Psi)
    }
    if (!is.null(data)) {
      data <- na.omit(data)
      cov_mat <- var(data, y = NULL, na.rm = TRUE)
      pop_lambda <- pop_psi <- NULL
    }
    if (!is.null(true_cov_matrix)) {
      if (!isSymmetric(true_cov_matrix, tol = 1e-05)) {
        .dmar_refuse_argument(
          "Input a symmetric covariance or correlation matrix 'true_cov_matrix'",
          arg = "true_cov_matrix")
      }
      .ss_aipe_check_S(true_cov_matrix)
      cov_mat <- true_cov_matrix
      n_items <- dim(true_cov_matrix)[1]
      pop_lambda <- pop_psi <- NULL
    }
    # With two items the single-factor model with free loadings is not
    # identified (three moments, four parameters), and its information
    # matrix is singular.
    if (ncol(cov_mat) < 3) {
      .dmar_refuse_argument(c(
        "The congeneric model needs at least three items: with two, the ",
        "single-factor model with free loadings is not identified. Use ",
        "the True Score model for two items."), arg = NULL)
    }
    pop <- .ss_aipe_population(cov_mat, pop_lambda, pop_psi)
    # A loading of zero leaves its item out of the factor, and with fewer
    # than three items on it the model is not identified at the population
    # (its information matrix is singular), whatever the number of columns.
    if (sum(pop$lambda != 0) < 3) {
      .dmar_refuse_argument(c(
        "The congeneric model needs at least three items with a nonzero ",
        "loading: with fewer, the single-factor model with free loadings ",
        "is not identified at the population, and omega has no sampling ",
        "variance to plan with."), arg = NULL)
    }
    # The delta method variance of the model implied omega from the
    # information of the single-factor model at the population, the
    # variance behind reliability_omega(denominator = "model_implied",
    # ci_method = "ml"); see .ss_aipe_omega_moments(). Before 1.1.0 this
    # was Terry and Kelley's (2012) Eq. 19, the variance of alpha, which
    # they used only as a starting value for the congeneric model (p. 390).
    mom <- .ss_aipe_omega_moments(cov_mat, pop$lambda, pop$psi)
    Crit_Value <- qnorm(1 - (1 - conf_level) / 2)
    Nec_N <- mom$V / (width / (2 * Crit_Value))^2 + 1
    # Items that load opposite to the others (.ss_aipe_opposite_items()).
    # Without an assurance the closed form is the plan, and it is refused
    # where its delta method is degenerate (.ss_aipe_check_delta()); an
    # assurance search counts the widths themselves and is not refused.
    opposite <- .ss_aipe_opposite_items(pop$lambda, cov_mat)
    if (is.null(assurance)) {
      .ss_aipe_check_delta(mom, width, Crit_Value, opposite, fmi)
    }
    if (length(opposite)) {
      .dmar_warn(.ss_aipe_opposite_message(opposite), "result_qualified")
    }
  }
  # The number of items the Monte Carlo draws need. A population given by
  # data never set it, and an assurance plan from data stopped in rep(0, n_items)
  # before 1.1.0.
  n_items <- ncol(cov_mat)
  feasibility_floor <- is.null(assurance) && Model_to_Use == "Congeneric" &&
    Nec_N < n_items + 1
  if (feasibility_floor) {
    Nec_N <- n_items + 1
    .dmar_warn(c(
      "The congeneric closed form is below items + 1 cases. The complete-data ",
      "plan is floored at ", n_items + 1, ", the minimum that permits a full-rank ",
      "sample covariance matrix. This is only a feasibility minimum, not ",
      "evidence that the width target is attained. Use an assurance plan ",
      "to assess the intervals at small sample sizes."), "result_qualified")
  }
  if (!is.null(assurance)) {
    if (assurance > 1) {
      assurance <- assurance / 100
    }
    initial_assurance_N <- ceiling(Nec_N) + 1
    # The search never goes below the smallest sample whose covariance
    # matrix can be nonsingular.
    n_min <- n_items + 1
    # Nor does it start there with the factor analytic interval, whose
    # single-factor fit needs a nonsingular sample covariance matrix: a
    # closed form below n_items cases (a highly reliable population planned for a
    # wide interval) would start the search at a size where no sample
    # returns an interval, and it stopped there. The normal theory
    # intervals need only alpha and the sample matrix, and start where they
    # always did.
    if (Type_to_Use == "Factor Analytic") {
      initial_assurance_N <- max(initial_assurance_N, n_min)
    }
    # Under missing data the search runs at the complete data (effective)
    # size, from the complete data closed form above, and the answer is
    # mapped back at the end. A starting size the user gives is a size
    # under missing data, so it is converted to its effective size first.
    if (!is.null(start_N) && fmi > 0) {
      start_N <- .aipe_ceiling(.aipe_fmi_effective(start_N, fmi))
    }
    # One Monte Carlo pass at the size n: n_iter samples of size n are drawn
    # from the population, and the interval of `type` is computed on each.
    # The assurance is the proportion of the n_iter samples whose interval
    # is no wider than `width` (Terry & Kelley, 2012, p. 373: the
    # probability of a confidence interval no wider than desired). A sample
    # for which the estimation functions return no interval
    # (.ss_aipe_reliability_ci() stops) did not give a narrow one, so it is
    # counted as wide, its slot holding Inf; through 1.1.0's development it
    # was set aside and the proportion taken over the rest, which with many
    # such samples (a third of them at 78 cases for loadings .8, .3, .2,
    # and .2) overstated the assurance. The draws stop as soon as the
    # outcome is settled: the slots not yet drawn hold `fill`, 0 (an
    # interval counted as narrow) when a pass is looking for a shortfall and
    # 1 (counted as wide) when it is looking for success, so the proportion
    # can be read before the last draw. With fill 0 the pass meets the
    # target when the proportion is at least `assurance`; with fill 1, when
    # it is above it. These are the rules of the search before 1.1.0, kept
    # so that a search that never steps down, and in which every sample
    # returns an interval, draws exactly what it drew then. Until some
    # sample has returned an interval the pass does not stop early, so a
    # pass in which none does is recognized: it reports that (`none`) with
    # the first error, rather than the missing value the earlier search
    # stopped on. `no_interval` counts the samples drawn that returned no
    # interval.
    mc_pass <- function(n, n_iter, fill, stage) {
      CI_Result <- rep(fill, n_iter)
      first_error <- NULL
      n_drawn <- 0L
      no_interval <- 0L
      Difference <- NA_real_
      for (iters in seq_len(n_iter)) {
        sim_data <- var(MASS::mvrnorm(n = n, mu = rep(
          0,
          n_items
        ), Sigma = cov_mat, tol = 1e-06, empirical = FALSE))
        CI_Result_raw <- try(.ss_aipe_reliability_ci(
          S = sim_data,
          N = n, model = Model_to_Use, type = Type_to_Use,
          conf_level = conf_level
        ), silent = TRUE)
        n_drawn <- iters
        if (inherits(CI_Result_raw, "try-error")) {
          CI_Result[iters] <- Inf
          no_interval <- no_interval + 1L
          if (is.null(first_error)) {
            first_error <- conditionMessage(attr(CI_Result_raw, "condition"))
          }
        } else {
          CI_Result[iters] <- CI_Result_raw$upper -
            CI_Result_raw$lower
        }
        Difference <- mean(CI_Result < width) - assurance
        if (verbose == TRUE) print(paste("The current assurance is", mean(CI_Result < width), "at the current sample size of", n, "; you are in the", sQuote(stage, FALSE), "stage."))
        if (no_interval == n_drawn) next
        if (fill == 0 && Difference < 0) break
        if (fill == 1 && Difference > 0) break
      }
      none <- no_interval == n_drawn
      list(difference = Difference, none = none,
           pass = !none && (if (fill == 0) Difference >= 0 else Difference > 0),
           empirical = mean(CI_Result < width),
           widths = CI_Result, first_error = first_error,
           no_interval = no_interval)
    }
    # A size at which no sample returned an interval cannot be judged, and
    # going on to the next size would not help.
    stop_none <- function(n, res) {
      .dmar_refuse_argument(c(
        "None of the Monte Carlo samples of size ", n, " returned a ",
        "confidence interval, so the assurance cannot be estimated. The ",
        "first failure: ", res$first_error), arg = NULL)
    }
    # The initial stage starts from the closed form, which targets the
    # expected width. When that size falls short, the size is raised one
    # case at a time until a pass meets the target, as before 1.1.0. When
    # it already meets the target (an assurance near or below one half, or
    # an interval narrower than the one the closed form describes, as the
    # factor analytic interval of the parallel model can be at a population
    # that is not parallel), the search
    # goes down rather than returning the starting size, which before 1.1.0
    # it did whatever the size the target needed. A pass that meets the
    # target runs every iteration, so going down one case at a time is
    # slow; the size is instead lowered to where the widths of the passing
    # sample would just meet the target, since the width of a Wald interval
    # shrinks as 1 / sqrt(N - 1): with q the `assurance` quantile of the
    # widths at the size n, the next size tried is (n - 1) (q / width)^2 + 1.
    # The step is repeated while the size tried meets the target. Once one
    # falls short, the size is raised one case at a time from there, so the
    # stage ends, as it does going up, at the first size that meets the
    # target above one that did not, and the final stage is the same.
    n_i <- initial_assurance_N
    if (is.null(start_N)) {
      res <- mc_pass(n_i, G_initial, 0, "initial")
      if (res$none) stop_none(n_i, res)
      if (res$pass) {
        while (n_i > n_min) {
          # A sample with no interval holds Inf; a pass that met the target
          # has at most 1 - assurance of them, so q is finite or, at the
          # boundary, Inf, which lowers the size by one case.
          q <- stats::quantile(res$widths, probs = assurance, names = FALSE)
          if (!is.finite(q)) q <- width
          n_try <- max(min(ceiling((n_i - 1) * (q / width)^2 + 1), n_i - 1),
                       n_min)
          res_try <- mc_pass(n_try, G_initial, 0, "initial")
          if (res_try$pass) {
            n_i <- n_try
            res <- res_try
            next
          }
          repeat {
            n_try <- n_try + 1
            if (n_try >= n_i) break
            res_try <- mc_pass(n_try, G_initial, 0, "initial")
            if (res_try$pass) {
              n_i <- n_try
              break
            }
          }
          break
        }
      } else {
        repeat {
          n_i <- n_i + 1
          res <- mc_pass(n_i, G_initial, 0, "initial")
          if (res$none) stop_none(n_i, res)
          if (res$pass) break
        }
      }
    }
    if (!is.null(start_N)) {
      n_i <- start_N
      res <- mc_pass(n_i, G, 0, "final")
      if (res$none) stop_none(n_i, res)
      if (res$difference > 0) {
        repeat {
          res <- mc_pass(n_i, G, 1, "final")
          if (!res$pass || n_i <= n_min) break
          n_i <- n_i - 1
        }
      }
    }
    repeat {
      res <- mc_pass(n_i, G, 0, "final")
      if (res$none) stop_none(n_i, res)
      if (res$pass) break
      n_i <- n_i + 1
    }
    Nec_N_assurance <- n_i
    assurance_attained <- res$empirical
    no_interval <- res$no_interval
  }
  # Nec_N (closed form) and Nec_N_assurance (search) are complete data
  # sizes. With fmi > 0 each is mapped to the size under missing data that
  # carries the same information, N = (N_c - 1) / (1 - fmi) + 1 (see
  # R/ss_aipe_missing_internals.R), and the complete data answer travels as
  # an attribute. With fmi = 0 nothing is mapped and no attribute is added,
  # so the result is the one the function returned before fmi existed.
  N_complete <- if (is.null(assurance)) Nec_N else Nec_N_assurance
  N_plan <- if (fmi > 0) {
    .aipe_ceiling(.aipe_fmi_inflate(N_complete, fmi))
  } else {
    ceiling(N_complete)
  }
  # Both paths carry the AIPE planner class, so tidy() and glance() give the
  # planned size whichever path answered (before 1.1.0 the assurance path
  # fell through to the generic table methods); the assurance table echoes
  # the width, and the closed form's single row does not, so its tidy()
  # width is NA, as in 1.0.1.
  out <- if (is.null(assurance)) {
    .as_dmar_tbl(data.frame(term = "necessary_N", value = N_plan), conf_level = conf_level, subclass = "dmar_ss_aipe")
  } else {
    .as_dmar_tbl(data.frame(
      term = c("necessary_N", "width", "assurance", "assurance_attained", "G"),
      value = c(N_plan, width, assurance, assurance_attained, G)
    ), conf_level = conf_level, subclass = "dmar_ss_aipe")
  }
  if (fmi > 0) {
    attr(out, "fmi") <- fmi
    attr(out, "necessary_N_complete") <- ceiling(N_complete)
  }
  if (feasibility_floor) attr(out, "feasibility_floor") <- n_items + 1L
  # The samples of the final pass that returned no interval, counted as
  # wide in assurance_attained, are recorded and named under the table
  # (.ss_aipe_no_interval_line()), so a plan at which some studies would
  # yield no coefficient says how many.
  if (!is.null(assurance)) {
    attr(out, "no_interval") <- as.integer(no_interval)
  }
  out
}

# A Parallel population has equal variances and equal off-diagonal
# covariances in the original item units. A common rescaling leaves this
# check unchanged; standardizing each item would change the model.
.ss_aipe_check_parallel <- function(S) {
  scale <- max(abs(S))
  tol <- sqrt(.Machine$double.eps) * scale
  d <- diag(S)
  off <- S[lower.tri(S)]
  if (diff(range(d)) > tol || diff(range(off)) > tol) {
    .dmar_refuse_argument(c(
      "The supplied population is not Parallel: item variances and ",
      "off-diagonal covariances must each be equal within a tolerance of ",
      "sqrt(machine precision) times the largest absolute covariance. ",
      "Use the True Score model for equal loadings with unequal error ",
      "variances, or the Congeneric model for unequal loadings. A sample ",
      "covariance from data need not satisfy an exact population model; ",
      "supply the intended parallel population through 'true_factor_loading' and ",
      "'true_residual_var' when that is the planning assumption."),
      arg = c("true_factor_loading", "true_residual_var"))
  }
  invisible(NULL)
}


# .ss_aipe_reliability_ci(S, N, model, type, conf_level)
#
# Confidence interval for the reliability coefficient implied by a sample
# covariance (or correlation) matrix, used by the Monte Carlo assurance
# search in ss_aipe_reliability(). The caller has already normalized
# `model` to one of "Parallel", "True Score", or "Congeneric" and `type`
# to "Factor Analytic" or "Normal Theory".
#
# The computation is the package's own reliability machinery, the same
# engine the estimation side of DMAR uses:
#
#   - "Factor Analytic": .omega_fit_cfa() fits the single-factor model by
#     maximum likelihood (equal loadings for the parallel and tau
#     equivalent models, free loadings for the congeneric model) and
#     .ci_omega_delta() builds the delta method Wald interval on the
#     omega defined parameter. This is the interval of
#     reliability_omega(denominator = "model_implied", ci_method = "ml")
#     and of reliability_alpha(estimator = "model_implied",
#     ci_method = "ml"). Note that the factor analytic treatment of the
#     parallel model matches the historical planner engine: equal
#     loadings are imposed but the error variances stay free, so the
#     parallel and tau equivalent models share the same factor analytic
#     interval.
#
#   - "Normal Theory": the van Zyl, Neudecker, and Nel (2000) closed
#     form standard error from .ci_alpha_ml() for the tau equivalent
#     model, and its compound symmetry simplification,
#     Var(alpha_hat) = 2 (1 - alpha)^2 J / ((N - 1) (J - 1)), for the
#     parallel model. The congeneric model has no normal theory form
#     and is an error.
#
# A sample is treated as the estimation functions treat it, so that the
# search counts the intervals a user would be given:
#
#   - A factor analytic fit that yields no coefficient stops, and the
#     planner counts the sample as one whose interval is wider than the
#     target (mc_pass() in ss_aipe_reliability()): a fit that did not
#     converge or stopped with an error, three items with no single-factor
#     solution, and a solution that has run away, whose error variances
#     sum below zero so that the model implied coefficient is above one
#     (.omega_fit_verdict()). reliability_omega() and the model implied
#     reliability_alpha() refuse the last two (.omega_point_gate()).
#     Before DMAR 1.1.0 the search kept a runaway solution, as the
#     estimation functions then did: with a standard error lavaan could
#     not compute (0) its interval came back inverted, an estimate of 104
#     with limits 104 and 1, and counted as narrow, and with a large one
#     it counted as wide. An improper solution whose error variances
#     still sum to a nonnegative number is kept, as those functions keep
#     it (with a warning).
#   - Each limit is clamped into [0, 1], as .relia_result() clamps the
#     limits the estimation functions return, so the limits never come
#     back inverted: a negative alpha whose normal theory interval lies
#     wholly below zero has the interval [0, 0], which is what
#     reliability_alpha() reports, rather than a lower limit of 0 above a
#     negative upper limit.
#
# Inputs:
#   S          : symmetric covariance or correlation matrix among items.
#   N          : sample size that generated S.
#   model      : "Parallel", "True Score", or "Congeneric".
#   type       : "Factor Analytic" or "Normal Theory".
#   conf_level : confidence level.
#
# Returns:
#   A list with `estimate`, `se`, `lower`, and `upper`. Stops (rather
#   than returning NA) when the factor analytic fit yields no
#   coefficient; the planner counts that sample as a wide interval.
.ss_aipe_reliability_ci <- function(S, N, model, type, conf_level) {
  clamp <- function(x) min(max(x, 0), 1)
  if (type == "Factor Analytic") {
    equal_loadings <- model %in% c("Parallel", "True Score")
    fit <- .omega_fit_cfa(S = S, N = as.numeric(N),
                          equal_loadings = equal_loadings, se = "standard")
    verdict <- .omega_fit_verdict(fit)
    # A sample whose likelihood has no maximum has run away along a ridge,
    # and reliability_omega() refuses it (.lavaan_no_maximum()).
    if (!is.null(fit$no_maximum)) {
      .dmar_refuse_argument(c(
        "The single-factor solution has run away (",
        .no_max_short(fit$no_maximum), "), and reliability_omega() ",
        "returns no coefficient for it."), arg = NULL)
    }
    if (verdict == "not_converged") {
      .dmar_refuse(
        if (is.null(fit$fit_error)) {
     "The single-factor model for coefficient omega did not converge."
   } else {
     paste0("The single-factor model for coefficient omega could not be ",
            "fit: ", fit$fit_error)
   }, "not_converged")
    }
    if (verdict == "no_solution") {
      .dmar_refuse_argument(c(
        "The three items have no single-factor solution (the product of ",
        "their covariances is not positive)."), arg = NULL)
    }
    if (verdict == "negative_error") {
      .dmar_refuse_argument(c(
        "The single-factor solution has run away: its error variances ",
        "sum below zero, so the model implied coefficient is above one, ",
        "and reliability_omega() returns no coefficient for it."), arg = NULL)
    }
    cl <- .ci_omega_delta(fit, conf_level = conf_level, logistic = FALSE)
    return(list(estimate = fit$omega, se = cl$se,
                lower = clamp(cl$lower), upper = clamp(cl$upper)))
  }
  if (type == "Normal Theory") {
    if (model == "Congeneric") {
      .dmar_refuse_argument(
        "The Congeneric model can not be used with `Normal Theory.'",
        arg = NULL)
    }
    p <- ncol(S)
    estimate <- .alpha_from_S(S)
    if (model == "Parallel") {
      crit <- stats::qnorm(1 - (1 - conf_level) / 2)
      se <- sqrt((2 * (1 - estimate)^2 * p) / ((N - 1) * (p - 1)))
      cl <- list(se = se,
                 lower = estimate - crit * se,
                 upper = estimate + crit * se)
    } else {
      cl <- .ci_alpha_ml(estimate, S = S, N = N, J = p,
                         conf_level = conf_level, logistic = FALSE)
    }
    return(list(estimate = estimate, se = cl$se,
                lower = clamp(cl$lower), upper = clamp(cl$upper)))
  }
  .dmar_refuse_argument(c("Unknown interval type: ", type), arg = NULL)
}


# .ss_aipe_omega_avar(S, lambda = NULL, psi = NULL, equal_loadings = FALSE)
#
# The asymptotic variance of the model implied coefficient omega, times the
# sample size, at a population covariance matrix S under the single-factor
# model: the closed form of ss_aipe_reliability() for the congeneric model
# (free loadings) and, with equal_loadings = TRUE, for the factor analytic
# interval of the tau equivalent model, whose expected interval width at N
# is taken to be 2 z sqrt(V / (N - 1)). lavaan's standard error has N where
# this has N - 1; the planner uses N - 1 as its alpha closed forms do,
# which adds one case and keeps the fmi mapping (N - 1 scaled by
# 1 / (1 - fmi)) the same for every model.
#
# With the factor variance fixed at 1 (lavaan's std.lv), the parameters are
# the J loadings lambda and the J error variances psi, and
#
#   omega = L^2 / (L^2 + P),   L = sum(lambda),  P = sum(psi),
#
# the model implied coefficient of reliability_omega(denominator =
# "model_implied") and of Terry and Kelley (2012, Eq. 13). Its gradient is
# 2 L P / (L^2 + P)^2 for every loading and -L^2 / (L^2 + P)^2 for every
# error variance (Terry & Kelley, 2012, Eqs. 14 and 15), and the expected
# information of the normal model for one observation is
#
#   I_ab = 1/2 tr( Sigma^-1 dSigma/d theta_a  Sigma^-1 dSigma/d theta_b ),
#
# with dSigma/d lambda_j = e_j lambda' + lambda e_j' and dSigma/d psi_j =
# e_j e_j'. The delta method then gives V = g' I^-1 g (Raykov, 2002; Terry
# & Kelley, 2012, Eq. 16), and V / N is the squared standard error lavaan
# reports for the omega defined parameter in a maximum likelihood fit, the
# one reliability_omega(denominator = "model_implied", ci_method = "ml")
# uses. V does not change when S is multiplied by a constant.
#
# With equal loadings (the tau equivalent model, whose model implied
# coefficient is the one reliability_alpha(estimator = "model_implied",
# ci_method = "ml") estimates) the J loadings are one parameter, lambda_j =
# lambda for every j. Writing the J + 1 parameters through the 2J x (J + 1)
# matrix C that sends lambda to every loading and each psi_j to itself,
# the information is C' I C and the gradient C' g: the loading's
# information is the sum of the loading block of I, and its gradient the
# sum of the loadings' gradients. At a population of parallel items (equal
# error variances too) V is the compound symmetry variance of van Zyl,
# Neudecker, and Nel (2000), 2 (1 - omega)^2 J / (J - 1).
#
# The loadings and error variances are those supplied when the population
# was given by them. When it was given by a covariance matrix (S or data),
# the single-factor model is fit to it by maximum likelihood first
# (.ss_aipe_population()). The work is done by .ss_aipe_omega_moments(),
# which also returns the variance of the sum of the loadings.
.ss_aipe_omega_avar <- function(S, lambda = NULL, psi = NULL,
                                equal_loadings = FALSE) {
  .ss_aipe_omega_moments(S, lambda, psi, equal_loadings = equal_loadings)$V
}


# .ss_aipe_omega_moments(S, lambda = NULL, psi = NULL,
#                        equal_loadings = FALSE)
#
# The computation behind .ss_aipe_omega_avar() (see the comment there). It
# returns a list with
#
#   V      : the asymptotic variance of the model implied coefficient,
#            times the sample size, g' I^-1 g;
#   V_L    : the same for the sum of the loadings L, d' I^-1 d with d the
#            gradient of L (1 for each loading and 0 for each error
#            variance; J for the one loading of the equal loadings model),
#            which .ss_aipe_check_delta() reads;
#   L      : the sum of the loadings;
#   lambda, psi : the loadings and error variances used.
.ss_aipe_omega_moments <- function(S, lambda = NULL, psi = NULL,
                                   equal_loadings = FALSE) {
  pop <- .ss_aipe_population(S, lambda, psi, equal_loadings = equal_loadings)
  lambda <- pop$lambda
  psi <- pop$psi
  J <- length(lambda)
  Sigma <- tcrossprod(lambda) + diag(psi, nrow = J)
  Sigma_inv <- solve(Sigma)
  # Sigma^-1 dSigma/d theta_a for the 2J parameters, loadings first.
  A <- vector("list", 2L * J)
  for (j in seq_len(J)) {
    e <- numeric(J)
    e[j] <- 1
    A[[j]] <- Sigma_inv %*% (tcrossprod(e, lambda) + tcrossprod(lambda, e))
    A[[J + j]] <- Sigma_inv %*% tcrossprod(e)
  }
  info <- matrix(0, 2L * J, 2L * J)
  for (a in seq_len(2L * J)) {
    for (b in a:(2L * J)) {
      # tr(A_a A_b) as the sum of the elementwise product of A_a and t(A_b).
      info[a, b] <- info[b, a] <- 0.5 * sum(A[[a]] * t(A[[b]]))
    }
  }
  L <- sum(lambda)
  P <- sum(psi)
  grad <- c(rep(2 * L * P, J), rep(-L^2, J)) / (L^2 + P)^2
  grad_L <- c(rep(1, J), rep(0, J))
  # The parameters in the units of the items: each loading divided by its
  # item's standard deviation and each error variance by its item's
  # variance (the equal loading by the root mean item variance). V, a
  # quadratic form in the gradient and the inverse information, does not
  # depend on the units of the parameters, but the conditioning of the
  # information does, and in these units it says whether the model is
  # identified at the population whatever the units of the items.
  sd_item <- sqrt(diag(Sigma))
  scale <- c(sd_item, sd_item^2)
  if (equal_loadings) {
    C <- rbind(cbind(1, matrix(0, J, J)), cbind(0, diag(J)))
    info <- crossprod(C, info %*% C)
    grad <- as.numeric(crossprod(C, grad))
    grad_L <- as.numeric(crossprod(C, grad_L))
    scale <- c(sqrt(mean(sd_item^2)), sd_item^2)
  }
  info <- info * tcrossprod(scale)
  grad <- grad * scale
  grad_L <- grad_L * scale
  # A loading that is zero, or indistinguishable from zero, on all but two
  # items leaves the information singular: the two loadings are known only
  # through their product. The tolerance refuses standardized loadings
  # below about .0001, a population the plan could only put in the
  # millions; a fit to a matrix whose items load zero leaves them near
  # 1e-8.
  x <- tryCatch(solve(info, cbind(grad, grad_L), tol = 1e-8),
                error = function(e) NULL)
  if (is.null(x)) {
    .dmar_refuse_argument(c(
      "The single-factor model is not identified at this population: ",
      "its information matrix is singular, as it is when fewer than three ",
      "items have a nonzero loading, so the coefficient has no sampling ",
      "variance to plan with."), arg = NULL)
  }
  list(V = sum(grad * x[, 1L]), V_L = sum(grad_L * x[, 2L]), L = L,
       lambda = lambda, psi = psi)
}


# .ss_aipe_population(S, lambda = NULL, psi = NULL, equal_loadings = FALSE)
#
# The loadings and error variances of the planning population under the
# single-factor model: those supplied, when the population was given by
# them, and otherwise the maximum likelihood fit of the model (free
# loadings, or equal loadings for the tau equivalent model) to the
# population covariance matrix S, at an N large enough that the fit is the
# population's. A fit that does not converge, is improper (a negative error
# variance), or has no solution (three items) stops with a message, since
# the variance would not describe a proper population.
#
# Returns a list with `lambda` and `psi`, numeric vectors of length J.
.ss_aipe_population <- function(S, lambda = NULL, psi = NULL,
                                equal_loadings = FALSE) {
  if (!is.null(lambda) && !is.null(psi)) {
    return(list(lambda = as.numeric(lambda), psi = as.numeric(psi)))
  }
  fit <- suppressWarnings(.omega_fit_cfa(S = S, N = 1e6,
                                         equal_loadings = equal_loadings,
                                         se = "none"))
  if (!identical(.omega_fit_verdict(fit), "proper")) {
    .dmar_refuse(c(
      "The single-factor model",
      if (equal_loadings) " with equal loadings" else "",
      " fit to the population covariance matrix is not a proper solution ",
      "(it did not converge, has a negative error variance, or has no ",
      "solution), so the ",
      if (equal_loadings) "factor analytic" else "congeneric",
      " closed form cannot be computed from it. A population to plan for ",
      "is one the model describes with positive error variances; check ",
      "the matrix."), "computation_failed")
  }
  list(lambda = fit$loadings, psi = fit$errors)
}


# .ss_aipe_check_psi(psi_square, lambda = NULL, normal_theory = FALSE)
#
# The error variances of a population given by its loadings (the true_item_covariance
# and lambda routes of the True Score and Congeneric models). Each must be
# a finite number. Where a single-factor model is fit (the congeneric
# model, and the tau equivalent model with the factor analytic interval)
# each must also be above zero, as the Parallel model has always required
# of its one: a negative error variance is not a variance, and a zero one
# makes the item free of error, the boundary at which the single-factor
# solution is improper. The same population given as S is refused when its
# single-factor fit is improper (.ss_aipe_population()), so the two routes
# agree. Through 1.0.1 both planned it without a word: loadings .7, .6, .5,
# and .5 with error variances .51, -.1, .75, and .75 gave 81 cases for the
# congeneric model (and, as S, for the tau equivalent one), and one loading
# of .7 on all four items with those error variances gave 46 for the tau
# equivalent model whatever the type.
#
# With the normal theory interval of the tau equivalent model
# (normal_theory = TRUE) nothing is fit: the loading and the error
# variances only build the covariance matrix the plan is made from, as S
# gives it directly, and the variance of alpha of van Zyl, Neudecker, and
# Nel (2000) needs only a covariance matrix. That matrix, lambda lambda' +
# diag(psi_square), must then be a covariance matrix, as S must
# (.ss_aipe_check_S()), and an error variance at or below zero that leaves
# it one is planned, as the same matrix given as S is: one loading of .7
# with error variances .51, -.1, .75, and .75 gives 74 cases by either
# route, and with 0, .5, .5, and .5, 110. The matrix is checked before
# covmat_from_cfa() builds it, which would warn first about its
# determinant.
.ss_aipe_check_psi <- function(psi_square, lambda = NULL,
                               normal_theory = FALSE) {
  if (!is.numeric(psi_square) || !all(is.finite(psi_square))) {
    .dmar_refuse_argument(
      "Every value of 'true_residual_var' must be a finite number.",
      arg = "true_residual_var")
  }
  if (normal_theory) {
    .ss_aipe_check_S(tcrossprod(lambda) + diag(psi_square,
                                                nrow = length(psi_square)),
                     implied = TRUE)
    return(invisible(NULL))
  }
  if (any(psi_square <= 0)) {
    .dmar_refuse_argument(c(
      "Every value of 'true_residual_var' must be greater than zero: an error ",
      "variance below zero is not a variance, and one of zero makes the ",
      "item free of error, a boundary at which the single-factor solution ",
      "is improper. Neither describes a population to plan for."),
      arg = "true_residual_var")
  }
}


# .ss_aipe_check_S(S, implied = FALSE)
#
# A population given as S must be a covariance matrix, positive
# semidefinite up to rounding: a negative eigenvalue means that some
# weighted sum of the items has a negative variance. The single-factor fit
# of the congeneric route refused such a matrix with a message about an
# improper solution; the alpha routes planned from it, and an assurance
# search stopped in MASS::mvrnorm() with "'Sigma' is not positive
# definite". With implied = TRUE the matrix is the one a loading and
# error variances imply on the normal theory route of the tau equivalent
# model (.ss_aipe_check_psi()), and the message names them.
.ss_aipe_check_S <- function(S, implied = FALSE) {
  ev <- .relia_negative_eigenvalue(S)
  if (!is.null(ev)) {
    .dmar_refuse_argument(c(
      if (implied) {
   "The matrix that 'true_factor_loading' (or 'true_item_covariance') and 'true_residual_var' imply"
 } else "'true_cov_matrix'", " is not a covariance matrix: its smallest eigenvalue ",
      if (isTRUE(attr(ev, "standardized"))) "after standardizing item units " else "after common rescaling ", "is ", .ss_aipe_num(ev), ", so some weighted ",
      "sum of the items would have a negative variance, which no ",
      "population has."), arg = NULL)
  }
}


# .ss_aipe_check_common_factor(cov_mat, lambda = NULL, psi = NULL,
#                              fitted = TRUE, parallel = FALSE)
#
# The factor analytic interval of the parallel and tau equivalent models is
# that of the coefficient implied by the single-factor model with equal
# loadings, lambda^2 J^2 / (lambda^2 J^2 + sum(psi)). Items to which that
# model gives a loading of zero share no common factor in it: the
# coefficient is zero, on the boundary of its range, and the closed form
# does not describe the interval. Through 1.0.1 such a population was
# planned (S with -.2 off the diagonal: 36880 cases at a width of .10 for
# the parallel model, and 4099 at .30 for both models; true_item_covariance = 0 gave
# 4099 at .10).
#
# For a population given by its loading (`fitted` FALSE; the true_item_covariance and
# lambda routes), the loading is the population's, and it is zero exactly
# when alpha of the matrix is zero or below. For one given as S or data
# the model is fit to it by maximum likelihood (.ss_aipe_population()),
# and the fitted loading decides, to the tolerance of
# .ss_aipe_opposite_items(), 1e-4 on standardized items. Alpha decides
# only for items of equal variance: at a loading of zero, with each error
# variance at its item's variance, the derivative of the ML discrepancy
# with respect to lambda^2 is minus the sum over pairs of s_ij / (s_ii
# s_jj), a sum weighted by the inverse variances, while the sign of alpha
# is that of the unweighted sum. Through the first versions of 1.1.0 alpha
# decided for every route. Three items of variance 1 that correlate .3
# and a fourth of variance 4 that correlates -.25 with them have alpha
# -0.276 but a fitted loading of .440 (coefficient .315, the interval
# reliability_alpha(estimator = "model_implied", ci_method = "ml") gives
# for them), and were refused; two items of variance 10 that covary 4 and
# two of variance 1 that covary -.3 have alpha .336 but a fitted loading
# of zero, and the tau equivalent model stopped with the message for an
# unidentified model while the parallel model planned 202 cases.
#
# The parallel model's closed form is the compound symmetry form in alpha
# of the population matrix (van Zyl, Neudecker, & Nel, 2000), which
# describes the
# factor analytic interval only for a population near the parallel one; with
# alpha at or below zero it is refused (`parallel`), and when the fitted
# coefficient is positive the message says the population is far from
# parallel and names the tau equivalent model, whose closed form is the
# fitted coefficient's. A fit that is improper leaves the parallel model
# to alpha, as before; the tau equivalent model stops on it in
# .ss_aipe_population().
#
# A coefficient above zero, however small, is planned: the loading enters
# the covariance matrix only through its square, of which the coefficient
# is a smooth function, so its variance does not vanish, and near zero the
# interval the estimation functions report, clamped at zero, only makes
# the plan conservative (four parallel items that covary .03 at a width of
# .30: 362 cases, at which the interval averages .236 wide).
#
# Returns the population's loadings and error variances (list(lambda,
# psi)), invisibly, or NULL when they were not computed; stops otherwise.
.ss_aipe_check_common_factor <- function(cov_mat, lambda = NULL, psi = NULL,
                                         fitted = TRUE, parallel = FALSE) {
  J <- ncol(cov_mat)
  alpha <- (J / (J - 1)) * (1 - sum(diag(cov_mat)) / sum(cov_mat))
  advice <- paste0(
    "Check the matrix (an item keyed in the opposite direction to the ",
    "others and not reverse scored makes its covariances negative), or use ",
    "type = \"Normal Theory\" to plan for alpha itself, whose estimate can ",
    "fall below zero.")
  boundary <- paste0(
    "the coefficient its factor analytic interval estimates is zero, on ",
    "the boundary of its range, where the closed form does not describe ",
    "that interval.")
  if (!fitted) {
    if (alpha > 0) {
      return(invisible(if (!is.null(lambda)) list(lambda = lambda, psi = psi)))
    }
    .dmar_refuse_argument(c(
      "The items of this population covary zero or less on average ",
      "(alpha of the population matrix is ", .ss_aipe_num(alpha), "), so ",
      "they share no common factor: the single-factor model with equal ",
      "loadings gives them a loading of zero, and ", boundary, " ",
      advice), arg = NULL)
  }
  pop <- if (parallel) {
    tryCatch(.ss_aipe_population(cov_mat, equal_loadings = TRUE),
             error = function(e) NULL)
  } else {
    .ss_aipe_population(cov_mat, lambda, psi, equal_loadings = TRUE)
  }
  if (!is.null(pop)) {
    std <- abs(pop$lambda[1L]) / sqrt(mean(diag(cov_mat)))
    if (std < 1e-4) {
      .dmar_refuse_argument(c(
        if (alpha <= 0) {
     paste0("The items of this population covary zero or less on ",
            "average (alpha of the population matrix is ",
            .ss_aipe_num(alpha), "), so they share no common factor: the ",
            "single-factor model with equal loadings gives them a ",
            "loading of zero, and ")
   } else {
     paste0("The single-factor model with equal loadings, fit to this ",
            "population by maximum likelihood, gives its items a loading ",
            "of zero, so they share no common factor in that model, ",
            "although alpha of the population matrix is ",
            .ss_aipe_num(alpha), ": the fit weighs each covariance by the ",
            "inverse of its two items' variances, alpha does not, and ",
            "with items of unequal variance the two can disagree. So ")
   }, boundary, " ", advice), arg = NULL)
    }
  }
  if (parallel && alpha <= 0) {
    omega <- if (is.null(pop)) NA_real_ else {
      L2 <- (J * pop$lambda[1L])^2
      L2 / (L2 + sum(pop$psi))
    }
    v <- range(diag(cov_mat))
    .dmar_refuse_argument(c(
      "Alpha of this population matrix is ", .ss_aipe_num(alpha),
      ", at or below zero, and the closed form of the parallel model is ",
      "a function of alpha, which it does not describe there.",
      if (is.finite(omega)) {
        paste0(" The single-factor model with equal loadings does give ",
               "the items a common factor (its coefficient, the one the ",
               "factor analytic interval estimates, is ",
               .ss_aipe_num(omega), "), but a population whose items ",
               "covary zero or less on the unweighted sum while that ",
               "coefficient is positive is far from parallel (its item ",
               "variances run from ", .ss_aipe_num(v[1L]), " to ",
               .ss_aipe_num(v[2L]), "). Plan with model = \"True ",
               "Score\", whose factor analytic closed form is that ",
               "interval's.")
      } else {
        paste0(" ", advice)
      }), arg = NULL)
  }
  invisible(pop)
}


# .ss_aipe_clamped_width_ratio(r, z)
#
# The mean width of the interval the estimation functions report for the
# model implied omega, over the width the delta method closed form plans
# for, when the sum of the loadings L is near zero relative to its
# standard error sigma_L at the planned N, r = sigma_L / |L|. Near zero
# omega_hat is about L_hat^2 / P (P the sum of the error variances) and its
# delta method standard error about 2 |L_hat| sigma_L / P. With X = L_hat /
# sigma_L normal with mean 1 / r and variance 1, the interval omega_hat
# -/+ z se, its lower limit clamped at zero as .relia_result() clamps it,
# has width (sigma_L^2 / P) min(4 z |X|, X^2 + 2 z |X|), and the closed form
# plans for 2 z times the population standard error, (sigma_L^2 / P) 4 z /
# r. The ratio is below one, the clamp making the plan conservative, for r
# up to 1.49 at 95 percent confidence (1.36 at 90 and 1.66 at 99 percent),
# and above one beyond: 1.21 at r = 2 and 1.86 at r = 3.35. Simulations of
# the factor analytic interval at the planned N (loadings .7, -.7, .5 and a
# fourth loading, 150 to 300 samples) gave .82, .82, .89, 1.21, and 2.2 at
# r = .41, 1.07, 1.52, 2.00, and 3.35 (at the planned N, rounded up), so
# the approximation errs a little toward refusing. The integral runs over
# 10 standard deviations either side of the mean of X.
.ss_aipe_clamped_width_ratio <- function(r, z) {
  if (!is.finite(r)) return(Inf)
  mu <- 1 / r
  f <- function(x) pmin(4 * z * abs(x), x^2 + 2 * z * abs(x)) *
    stats::dnorm(x, mean = mu)
  stats::integrate(f, mu - 10, mu + 10)$value / (4 * z * mu)
}


# .ss_aipe_check_delta(mom, width, crit, opposite, fmi)
#
# The congeneric closed form, the delta method variance of the model
# implied omega, takes omega to be normal with its first-order variance,
# which needs the sum of the loadings L to be large relative to its
# standard error. Where the loadings of opposite sign nearly cancel in the
# unweighted sum, it is not: omega then behaves as the square of a normal
# variable near zero, the gradient of omega shrinks with L, and the closed
# form understates the width of the interval, until at L = 0 it gives
# N = 1. The check evaluates r = sigma_L / |L| at the N the closed form
# gives, rounded up as the plan is, or at i + 1 cases when it gives fewer,
# the fewest at which the factor analytic interval exists (with i cases
# or fewer the sample covariance matrix of i items is singular, and the
# search of an assurance plan starts there too): r = sqrt(V_L / (N - 1)) /
# |L|. It refuses the plan when the reported interval would average wider
# than the target there (.ss_aipe_clamped_width_ratio() above one).
# Through the first versions of 1.1.0 r was taken at the unrounded
# closed form, whose N - 1 is a fraction of one case for a highly reliable
# population and a wide target (0.026 for four loadings of .99 and a
# width of .2), so r was large (4.4) for a reason that had nothing to do
# with L near zero, and the plan was refused with a message about omega's
# distribution near zero that does not apply near one, while 1.0.1
# returned 2 (four loadings of .95 at the same width were not refused).
# The ratio is derived for omega near zero, where L is small relative to
# its standard error; at a realizable N a population of one sign far from
# zero is not near that regime, and the check passes. A closed form below
# i + 1 is returned as it is, as in 1.0.x. Before this check the closed
# form planned 49 cases for loadings .7, -.7, .5, -.4 (omega .0038, r =
# 3.35), at which the interval averages .22 wide for a target of .10 and
# where about 95 are needed, and 1 for .7, -.7, .5, -.5 (1.0.x, with
# alpha's closed form: 12606 and 14078). With an assurance the closed form
# is only the start of a search that counts the widths themselves, so the
# check is not made. `opposite` (from .ss_aipe_opposite_items()) adds the
# advice to reverse score; `fmi` > 0 says the N is the complete data one.
.ss_aipe_check_delta <- function(mom, width, crit, opposite, fmi = 0) {
  J <- length(mom$lambda)
  N_cf <- ceiling(mom$V / (width / (2 * crit))^2 + 1)
  N_eval <- max(N_cf, J + 1)
  r <- if (mom$L == 0 || mom$V <= 0) {
    Inf
  } else {
    sqrt(mom$V_L / (N_eval - 1)) / abs(mom$L)
  }
  ratio <- .ss_aipe_clamped_width_ratio(r, crit)
  if (ratio <= 1) return(invisible(r))
  .dmar_refuse_argument(c(
    "The closed form cannot plan this interval. ",
    if (is.finite(r)) {
      paste0(if (N_eval > N_cf) {
        paste0("At ", N_eval, " cases, the fewest at which the interval ",
               "can be computed (the closed form gives ", N_cf, ")")
      } else {
        paste0("At the N it gives, ", N_cf, if (N_cf == 1) " case" else
                 " cases")
      }, if (fmi > 0) " with complete data" else "",
             ", the sum of the loadings (", .ss_aipe_num(mom$L),
             ") would be estimated with a standard error ",
             .ss_aipe_num(r), " times its size, so omega, which grows ",
             "with its square, is far from the normal distribution the ",
             "delta method takes, and the interval would average wider ",
             "than ", .ss_aipe_num(width), " (about ",
             .ss_aipe_num(ratio), " times as wide by a normal ",
             "approximation). ")
    } else {
      paste0("The sum of the loadings is zero, so omega is zero, its ",
             "gradient vanishes, and the delta method gives no sampling ",
             "variance (the closed form gives ", N_cf, " case",
             if (N_cf == 1) "" else "s", "). ")
    },
    if (length(opposite)) {
      paste0("The loadings differ in sign and cancel, or nearly so, in ",
             "the unweighted sum, as items keyed in opposite directions do ",
             "when none is reverse scored; reverse score such items before ",
             "planning, as before the coefficient is estimated. ")
    } else "",
    "To plan for this coefficient as it is, give an assurance: the ",
    "Monte Carlo search counts the widths of the intervals themselves ",
    "(an assurance of .5 plans for the median width)."), arg = NULL)
}


# .ss_aipe_opposite_items(lambda, Sigma)
#
# The items of a congeneric population that load in the direction opposite
# to the others: with the factor oriented so that the sum of the loadings
# is positive (by the largest loading when the sum is zero), the items whose
# standardized loading is below -1e-4, when some other item's is above
# 1e-4 (the tolerance keeps the loadings a fit leaves near zero, about
# 1e-8, out of it, as the identification tolerance of
# .ss_aipe_omega_moments() does). Returns their labels, the column names of
# Sigma or the items' positions, or character(0).
.ss_aipe_opposite_items <- function(lambda, Sigma) {
  std <- lambda / sqrt(diag(Sigma))
  orient <- sign(sum(lambda))
  if (orient == 0) orient <- sign(lambda[which.max(abs(lambda))])
  std <- orient * std
  opp <- which(std < -1e-4)
  if (!length(opp) || !any(std > 1e-4)) return(character(0))
  labels <- colnames(Sigma)
  if (is.null(labels)) labels <- as.character(seq_along(lambda))
  labels[opp]
}


# .ss_aipe_opposite_message(opposite)
#
# The warning ss_aipe_reliability() gives when items of a congeneric
# population load in the direction opposite to the others.
.ss_aipe_opposite_message <- function(opposite) {
  n <- length(opposite)
  items <- if (n == 1) {
    paste("item", opposite)
  } else if (n == 2) {
    paste("items", paste(opposite, collapse = " and "))
  } else {
    paste0("items ", paste(opposite[-n], collapse = ", "), ", and ",
           opposite[n])
  }
  paste0("The loadings of the population differ in sign: ", items,
         if (n == 1) " loads" else " load",
         " in the direction opposite to the other items, and so lower",
         if (n == 1) "s" else "", " omega of the unweighted sum, as an item ",
         "keyed in the opposite direction to the others does when it is not ",
         "reverse scored. The plan is for omega of the sum of the items as ",
         "they are given; reverse score such items before planning, as ",
         "before the coefficient is estimated.")
}


# .ss_aipe_num(x)
#
# A number for a message: three significant digits, without the padding
# formatC(format = "fg") leaves.
.ss_aipe_num <- function(x) trimws(formatC(x, digits = 3, format = "fg"))


# .ss_aipe_no_interval_line(x)
#
# The line written under an assurance plan of ss_aipe_reliability() whose
# final Monte Carlo pass had samples that returned no interval (the
# `no_interval` attribute), or NULL. Its caller is
# .dmar_tbl_footer.default() in R/dmar_output_helpers.R.
.ss_aipe_no_interval_line <- function(x) {
  k <- attr(x, "no_interval", exact = TRUE)
  if (!is.numeric(k) || length(k) != 1L || k <= 0) return(NULL)
  iter <- x$value[x$term == "G"]
  paste0("Of the ", format(iter, scientific = FALSE), " samples of the ",
         "final Monte Carlo pass, ", format(k, scientific = FALSE),
         " returned no interval (the single-factor fit did not converge, ",
         "had no solution, or had run away) and ", if (k == 1) "was" else
           "were", " counted as wider than the target width.")
}
