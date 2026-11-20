#' Information Criteria for Structural Equation Models
#'
#' Computes the information criteria used to choose among structural
#' equation models fitted to the same data: AIC (Akaike, 1974), BIC
#' (Schwarz, 1978), the sample size adjusted BIC (SABIC; Sclove, 1987),
#' Haughton's (1988) BIC (HBIC), and the information matrix based BIC
#' (IBIC) and scaled unit information prior BIC (SPBIC) that Bollen, Ray,
#' Zavisca, and Harden (2012) proposed and Bollen, Harden, Ray, and Zavisca
#' (2014) evaluated for structural equation models. Each criterion is minus
#' twice the maximized log likelihood plus a penalty for the free
#' parameters, smaller values are preferred, and only the difference
#' between two models fitted to the same data under the same likelihood
#' carries meaning. The function takes either a fitted \pkg{lavaan} model,
#' whose log likelihood, parameter count, sample size, and (for IBIC and
#' SPBIC) information matrices it reads, or those quantities as summary
#' input, so that a published comparison can be reproduced.
#' \code{\link{compare_sem}} compares a set of fitted models by one
#' criterion. A fitted model requires \pkg{lavaan} (in Suggests).
#'
#' @param fit A fitted \pkg{lavaan} model: a maximum likelihood fit of
#'   continuous variables under the normal likelihood (estimator ML or its
#'   robust variants MLM, MLMV, and MLR) that converged to an admissible
#'   solution, as \code{\link{sem_fit_measures}} judges it. See Details
#'   for what is refused. Supply either \code{fit} or the summary input
#'   below, not both.
#' @param log_likelihood For summary input, the maximized log likelihood
#'   \eqn{\ell}, including its constants (lavaan's \code{logl}).
#' @param npar For summary input, the number \emph{d} of independent free
#'   parameters of the likelihood whose maximum is \code{log_likelihood}, a
#'   positive whole number, counted the same way for every model to be
#'   compared. A fitted model supplies lavaan's count,
#'   \code{attr(logLik(fit), "df")}, which leaves out the saturated means of
#'   a model without a mean structure and, under \code{fixed.x = TRUE}, the
#'   moments of the exogenous observed variables; a difference between two
#'   models with the same settings is unaffected.
#' @param N For summary input, the sample size \emph{N} in the penalties of
#'   the BIC family, greater than 1. A fitted model supplies its total over
#'   groups.
#' @param measures The criteria to compute: any of \code{"AIC"},
#'   \code{"BIC"}, \code{"SABIC"}, \code{"HBIC"}, \code{"IBIC"}, and
#'   \code{"SPBIC"}, in either case. The default computes the four that need
#'   only \eqn{\ell}, \emph{d}, and \emph{N}. Every row of the table is
#'   reported whatever is requested, with \code{NA} for a criterion not
#'   computed.
#' @param parameter_estimate For SPBIC from summary input, the vector
#'   \eqn{\hat{\theta}} of the \emph{d} free parameter estimates, in the
#'   order of the rows of \code{information_observed}.
#' @param prior_mean For SPBIC, the vector \eqn{\theta^*} of prior means,
#'   in the order of \code{parameter_estimate}. Supply this or
#'   \code{zero_parameters}. With a fitted model and neither, the rule of
#'   Bollen et al. (2014, p. 6) applies: free loadings and regression
#'   coefficients are centered at zero and every other free parameter at its
#'   estimate (see Details).
#' @param zero_parameters For SPBIC, the positions or unique names of the
#'   entries of \code{parameter_estimate} to center at zero, every other
#'   entry being centered at its estimate. A saturated model needs an
#'   explicit choice; Bollen et al. (2014, p. 6) center its covariances at
#'   zero and its variances (and means) at their estimates.
#' @param information_observed For SPBIC from summary input, the observed
#'   information matrix \eqn{I_O(\hat{\theta})}: the negative Hessian of
#'   the log likelihood at the estimates, in the free parameters, in the
#'   order of \code{parameter_estimate}.
#' @param information_expected For IBIC from summary input, the expected
#'   (Fisher) information matrix \eqn{I_E(\hat{\theta})} at the estimates,
#'   in the free parameters. The inverse of a robust (sandwich) covariance
#'   matrix of the estimates is not an information matrix and is no
#'   substitute.
#' @param information_scale For summary input, whether the supplied
#'   information matrices are for the whole sample (\code{"total"}, the
#'   default; the inverse of the covariance matrix of the estimates under ML)
#'   or per observation (\code{"per_observation"}, what
#'   \code{lavaan::lavInspect()} reports), in which case they are multiplied
#'   by \code{N}. A fitted model supplies per observation matrices, which
#'   are multiplied by \emph{N}, and the argument is refused with one.
#' @param spbic_dimension Optional: the dimension used in the penalty of
#'   SPBIC, a positive whole number no greater than \code{npar}. By default
#'   it is \code{npar}, as in the general expression. Bollen et al.'s (2014)
#'   Table 6 used smaller counts for SPBIC, printed in parentheses (15 for
#'   the saturated model rather than 20), and states no rule for them; supply
#'   such a count only to reproduce that table, never to tune a ranking. It
#'   changes neither \code{parameter_estimate} nor the information
#'   matrices.
#'
#' @details
#' \subsection{The Criteria}{
#' Write \eqn{\ell} for the maximized log likelihood, \emph{d} for the
#' number of free parameters, and \emph{N} for the sample size. Each
#' criterion is \eqn{-2\ell} plus a penalty: \eqn{AIC = -2\ell + 2d}
#' (Akaike, 1974), \eqn{BIC = -2\ell + d\log N} (Schwarz, 1978),
#' \eqn{SABIC = -2\ell + d\log\{(N + 2)/24\}} (Sclove, 1987; lavaan's
#' \code{bic2}), and \eqn{HBIC = -2\ell + d\log\{N/(2\pi)\}} (Haughton,
#' 1988).
#'
#' Bollen et al. (2012) derive the BIC family from the Laplace
#' approximation to the log of a model's marginal likelihood,
#' \deqn{\log P(Y|M) = \ell + \log P(\hat{\theta}|M) + \frac{d}{2}\log(2\pi)
#'   - \frac{d}{2}\log N - \frac{1}{2}\log|\bar{I}_E(\hat{\theta})| +
#'   O(N^{-1/2})}{log P(Y|M) = l + log P(theta|M) + (d/2) log(2 pi) - (d/2)
#'   log N - (1/2) log|I_E(theta)/N| + O(N^(-1/2))}
#' (Bollen et al., 2014, Eq. 9), where \eqn{P(\hat{\theta}|M)} is the prior
#' density at the estimates and \eqn{\bar{I}_E = I_E/N} the expected
#' information per observation. Multiplied by \eqn{-2}, BIC keeps the first
#' and fourth terms, HBIC the third as well, and IBIC every term but the
#' prior: \eqn{IBIC = -2\ell + d\log\{N/(2\pi)\} +
#' \log|\bar{I}_E(\hat{\theta})|} (Bollen et al., 2012, Eq. 20; 2014, Eq.
#' 16), which equals \eqn{-2\ell - d\log(2\pi) + \log|I_E(\hat{\theta})|}
#' for the whole sample information. Recomputing the worked example of
#' Bollen et al. (2012, Table 1) from its printed inputs shows that the log
#' determinant it adds to \eqn{d\log\{N/(2\pi)\}} is that of the whole
#' sample information (its elements are about \eqn{N} times the per
#' observation ones the log likelihood implies), which counts \eqn{\log N}
#' twice relative to their Eq. 20; this function follows Eq. 20, which
#' reproduces the IBIC row of Bollen et al. (2014, Table 6).
#'
#' SPBIC replaces BIC's unit information prior with a normal prior centered
#' at \eqn{\theta^*} (\code{prior_mean}) whose variance, a multiple of the
#' inverse unit information, is chosen to maximize the marginal
#' likelihood, an empirical Bayes prior. With
#' \eqn{Q = (\hat{\theta} - \theta^*)' I_O(\hat{\theta})
#' (\hat{\theta} - \theta^*)} for the estimates \eqn{\hat{\theta}}
#' (\code{parameter_estimate}) and the whole sample observed information
#' \eqn{I_O}, \eqn{SPBIC = -2\ell + d\{1 + \log(Q/d)\}} when \eqn{Q > d}
#' (Bollen et al., 2012, Eq. 14; 2014, Eq. 18). When \eqn{Q \le d} the
#' maximizing prior variance is zero, the prior is a point mass at
#' \eqn{\theta^*}, and \eqn{SPBIC = -2\ell + Q} (Bollen et al., 2012, Eq.
#' 15); the two branches meet at \eqn{Q = d}. Bollen et al. (2014, note 7)
#' never met that branch and caution that its behavior is little known, and
#' the function warns when it applies.
#'
#' SABIC replaces \emph{N} in BIC's penalty with \eqn{(N + 2)/24}, the
#' penalty of Rissanen's (1978) shortest data description criterion as
#' Sclove (1987, criteria 4 and 5) gives it, which comes from the length
#' of a coded description of the data rather than from the Laplace
#' approximation; Lin, Huang, and Weng (2017) write it ABIC, and Mplus
#' and lavaan's \code{bic2} call it the sample size adjusted BIC. Its
#' penalty per parameter is smaller than BIC's, and it is not
#' positive when \eqn{N \le 22}, where the function warns. Lin et al. (2017)
#' found that it improved on BIC in selecting path models mainly with small
#' samples or small misspecifications, and that HBIC and SPBIC did best
#' overall. AIC estimates the relative expected Kullback-Leibler
#' discrepancy of a model rather than its marginal likelihood (Burnham &
#' Anderson, 2004). AICc is not offered: its correction is derived for the
#' univariate linear model and is not justified for structural equation
#' models in general.
#' }
#'
#' \subsection{Units and Identification}{
#' AIC, BIC, SABIC, and HBIC depend on a model only through \eqn{\ell} and
#' \emph{d}, which do not change when an observed variable is rescaled or a
#' latent variable is identified another way. IBIC and SPBIC also depend on
#' the information matrix in the free parameters as they are parameterized.
#'
#' IBIC keeps \eqn{\log|\bar{I}_E|} but drops the prior density of the
#' Laplace approximation, and that is the term whose Jacobian would cancel a
#' change of parameterization. When the units of a free parameter are
#' multiplied by \eqn{c^k}, \eqn{\log|\bar{I}_E|} changes by \eqn{-2k\log
#' c}. Rescaling an observed variable by \emph{c} therefore changes the
#' IBIC of one model relative to another by \eqn{-2k\log c} for each free
#' parameter the first has and the second lacks whose units carry that
#' variable's (\eqn{k = 1} for a loading or for a covariance with an
#' unscaled variable, \eqn{k = 2} for its variance), and the ranking of two
#' models fitted to the same data can reverse when one variable's units
#' change. Identifying the latent
#' variables by fixing their variances (\code{std.lv = TRUE}) rather than by
#' a marker indicator changes IBIC differences too, since the map between
#' the two parameterizations is nonlinear and its Jacobian differs between
#' models. An IBIC comparison is defined only for a stated choice of units
#' and identification. Standardizing every observed variable before fitting
#' every model, or keeping their original metrics, fixes the units, but it
#' is a convention rather than an invariant answer: report it, and read
#' IBIC beside BIC, HBIC, and SPBIC.
#'
#' SPBIC's \emph{Q} is invariant to a linear reparameterization when the
#' prior mean is transformed with the estimates, so rescaling an observed
#' variable does not change it under the default prior. That prior,
#' however, centers the free loadings and regression coefficients at zero
#' and the other free parameters at their estimates, and which parameters
#' are free loadings depends on the identification: under \code{std.lv =
#' TRUE} every loading is centered at zero, while under a marker indicator,
#' the convention of the example of Bollen et al. (2014), the scale sits in
#' the factor variance, which stays at its estimate. A mean structure adds
#' the saturated means to \emph{d}, and their information to
#' \eqn{\bar{I}_E}. Models compared by IBIC or SPBIC must therefore share
#' their identification and mean structure; \code{\link{compare_sem}}
#' requires every latent variable of every model to be identified the same
#' way (each by a marker indicator, or each by a fixed variance) and the
#' same mean structure, and refuses a set that does not, and for a fitted
#' model the \code{sem_information} attribute records both.
#' }
#'
#' \subsection{What a Fitted Model Supplies}{
#' \eqn{\ell} is lavaan's \code{logLik(fit)}, \emph{d} its \code{df}
#' attribute, and \emph{N} the total sample size. In a multiple group
#' model \eqn{\ell} is the sum over the groups and \emph{N} the total, as in
#' lavaan. Under \code{fixed.x = TRUE} (the default of
#' \code{lavaan::sem()}) or \code{conditional.x = TRUE}, \eqn{\ell} is the
#' log likelihood conditional on the exogenous observed variables, and
#' criteria compare only models with the same exogenous variables. Under
#' FIML (\code{missing = "ml"}) \eqn{\ell} is the observed data log
#' likelihood. Under \code{missing = "two.stage"} or
#' \code{"robust.two.stage"} it is the complete data normal log likelihood
#' at the stage one EM estimates of the means and covariances, so that AIC
#' and BIC are the two-stage criteria of Lai (2021, Eqs. 32 and 35), up to a
#' constant shared by all models fitted to the same data; Lai shows that,
#' unlike the FIML criteria, they select in large samples the model
#' complete data would select when the values are missing at random, and
#' that in samples of the usual size they miss it more often than complete
#' data would. They are comparable only with other two-stage fits of the same
#' variables with the same auxiliary variables in the first stage, and they
#' are the \code{AIC} and \code{BIC} rows \code{\link{fa_confirmatory}} reports
#' after two-stage estimation. For IBIC
#' and SPBIC, the expected and observed information matrices are lavaan's
#' per observation \code{lavInspect(fit, "information.expected")} and
#' \code{lavInspect(fit, "information.observed")}, multiplied by \emph{N};
#' they are the matrices the criteria are defined on under lavaan's defaults
#' \code{h1.information = "structured"} (the expected information at the
#' fitted moments) and \code{observed.information = "hessian"} (the negative
#' Hessian), and other settings are refused.
#'
#' A two-stage fit can have cases observed on an auxiliary variable
#' (lavaan's \code{aux}) but on none of the model's variables. lavaan's
#' first stage reads them, so they move the EM estimates, but its sample
#' size (\code{lavaan::lavInspect(fit, "nobs")}) counts them as empty and
#' leaves them out, and it forms the log likelihood as that sample size
#' times the complete data log likelihood per case at the EM estimates.
#' Its log likelihood, and the \emph{N} of the penalties of BIC, SABIC,
#' and HBIC, then rest on fewer cases than the moments. Both are computed
#' here from the fit computed again with every case the first stage read:
#' in one group, and in several groups that no constraint holds together,
#' the estimates do not depend on the count, and the log likelihood of
#' each group is its value per case times the cases its first stage read,
#' with \emph{N} their total. A line under the table says so, and the
#' \code{sem_information} attribute records the counts. On six Holzinger
#' and Swineford tests with the whole block missing for 77 of 301 pupils,
#' at random given general information (the auxiliary), lavaan's \emph{N}
#' for the two-factor model is 224 and its log likelihood -4040.25, which
#' is 224/301 of the -5429.09 reported here with \emph{N} = 301 (BIC
#' 10966.61 here, 8183.32 from lavaan's values). In several groups held
#' together by a constraint, such as equal loadings or intercepts, the
#' second stage weights each group by its count, so when the groups' counts
#' grow in different ratios the estimates move as well, and the log
#' likelihood is the one at the estimates with every case read (those of
#' semTools' \code{twostage()}), which the line says. On the same tests by
#' school with equal loadings and intercepts (130 and 94 cases with an
#' observed test, 156 and 145 read), it is -5393.34 with BIC 10957.89,
#' against -5393.52 and 10958.25 at lavaan's estimates and -4020.79 and
#' 8203.92 from lavaan's values.
#' \code{\link{fa_confirmatory}(missing = "two_stage")} and
#' \code{\link{measurement_invariance}} return fits that count the cases
#' already, and give the same criteria. FIML (\code{missing = "ml"}) keeps
#' the cases in its likelihood and in \emph{N} as well.
#' }
#'
#' \subsection{What Is Refused}{
#' A fitted model is refused when it did not converge or has an improper
#' solution, when it has ordered categorical indicators or was fitted by an
#' estimator other than ML (its fit function is then not a log likelihood),
#' when it uses the Wishart likelihood, and when it was fitted with sampling
#' weights: lavaan then maximizes a weighted pseudo log likelihood
#' (Asparouhov, 2005), for which neither the penalties nor the information
#' matrices have their likelihood justification. For a two-level model the
#' sample size in the BIC family's penalty is ambiguous (Bollen et al., 2012,
#' p. 295), so only AIC is computed from the fit, and the others need
#' summary input with the \emph{N} the user chooses. Automatic IBIC and
#' SPBIC further need a single group, single level ML fit (not MLM or MLR)
#' of complete data whose free parameters are unconstrained, with no
#' exploratory (rotated) blocks (their information matrices are singular)
#' and with the exogenous observed variables modeled (\code{fixed.x =
#' FALSE}). Each refusal says why and what to do; for another
#' parameterization, supply the quantities as summary input.
#' }
#'
#' \subsection{Relative to the Saturated Model}{
#' Bollen et al. (2014) report each criterion as the difference between the
#' hypothesized model and the saturated model, so that BIC becomes
#' \eqn{T - df\log N} for the likelihood ratio statistic \eqn{T} on
#' \eqn{df} degrees of freedom. Subtracting the saturated model's criterion
#' from the hypothesized model's gives that scale. A single model's
#' criterion cannot be reconstructed from \eqn{T} and \eqn{df} alone.
#' }
#'
#' @return A \code{dmar_tbl} (class \code{dmar_sem_information}) with
#'   columns \code{term} and \code{value} and, always in this order, the
#'   rows \code{AIC}, \code{BIC}, \code{SABIC}, \code{HBIC}, \code{IBIC},
#'   and \code{SPBIC} (\code{NA} for a criterion not requested),
#'   \code{log_likelihood}, \code{npar}, and \code{N}. The
#'   \code{sem_information} attribute records the source; for a fitted
#'   model the likelihood, the missing data method, the exogenous observed
#'   variables it conditions on, the mean structure, the number of groups,
#'   and, when IBIC or SPBIC is computed, the identification of each latent
#'   variable; the scale and provenance of the information matrices; the
#'   log determinant of the whole sample expected information for IBIC; for
#'   SPBIC the estimates, the prior means and the rule that set them,
#'   \emph{Q}, the penalty dimension and its rule, the penalty, and whether
#'   the boundary branch applied; for a two-stage fit whose first stage read
#'   cases observed on an auxiliary variable alone, \code{two_stage_rows},
#'   the cases read and \pkg{lavaan}'s count of them in each group, how many
#'   such cases there are, whether the fit was computed again with all of
#'   them (\code{corrected}), and whether its estimates depend on the
#'   counts, so that they were computed again too (\code{reestimated});
#'   and the lines of the printed footer.
#'
#' @references
#' Akaike, H. (1974). A new look at the statistical model identification.
#'   \emph{IEEE Transactions on Automatic Control, 19}(6), 716--723.
#'   \doi{10.1109/TAC.1974.1100705}
#'
#' Asparouhov, T. (2005). Sampling weights in latent variable modeling.
#'   \emph{Structural Equation Modeling, 12}(3), 411--434.
#'
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
#' Haughton, D. M. A. (1988). On the choice of a model to fit data from an
#'   exponential family. \emph{The Annals of Statistics, 16}(1), 342--355.
#'   \doi{10.1214/aos/1176350709}
#'
#' Lai, K. (2021). Using information criteria under missing data: Full
#'   information maximum likelihood versus two-stage estimation.
#'   \emph{Structural Equation Modeling, 28}(2), 278--291.
#'   \doi{10.1080/10705511.2020.1780925}
#'
#' Lin, L.-C., Huang, P.-H., & Weng, L.-J. (2017). Selecting path models in
#'   SEM: A comparison of model selection criteria. \emph{Structural
#'   Equation Modeling, 24}(6), 855--869.
#'   \doi{10.1080/10705511.2017.1363652}
#'
#' Rissanen, J. (1978). Modeling by shortest data description.
#'   \emph{Automatica, 14}(5), 465--471.
#'   \doi{10.1016/0005-1098(78)90005-5}
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. \emph{The Annals
#'   of Statistics, 6}(2), 461--464. \doi{10.1214/aos/1176344136}
#'
#' Sclove, S. L. (1987). Application of model-selection criteria to some
#'   problems in multivariate analysis. \emph{Psychometrika, 52}(3),
#'   333--343. \doi{10.1007/BF02294360}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{compare_sem}} to compare fitted models by one of
#'   these criteria, \code{\link{sem_fit_measures}} for the fit measures of
#'   one model, and the vignette
#'   \code{vignette("sem_fit_measures", package = "DMAR")}.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' # The three-factor model of nine of the holzinger_swineford
#' # tests, identified by a marker indicator (lavaan's default).
#' fit <- lavaan::cfa(
#'   "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
#'    textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
#'    speed   =~ t10_addition + t12_counting_groups_of_dots +
#'               t13_straight_and_curved_capitals",
#'   data = holzinger_swineford)
#' sem_information_criteria(fit)
#' sem_information_criteria(fit, measures = c("BIC", "IBIC", "SPBIC"))
#'
#' # Summary input: the autoregressive model of Bollen et al. (2014,
#' # Table 6), whose BIC relative to the saturated model is T - df log N.
#' ar <- sem_information_criteria(log_likelihood = -534.27 / 2, npar = 14,
#'                                N = 3912)
#' saturated <- sem_information_criteria(log_likelihood = 0, npar = 20,
#'                                       N = 3912)
#' ar$value[ar$term == "BIC"] - saturated$value[saturated$term == "BIC"]
#'
#' # SPBIC from summary input, with an explicit prior mean.
#' sem_information_criteria(log_likelihood = -120, npar = 2, N = 200,
#'                          measures = "SPBIC",
#'                          parameter_estimate = c(0.4, 1),
#'                          prior_mean = c(0, 1),
#'                          information_observed = diag(c(100, 50)))
#'
#' @export
sem_information_criteria <- function(fit = NULL, log_likelihood = NULL,
                                     npar = NULL, N = NULL,
                                     measures = c("AIC", "BIC", "SABIC", "HBIC"),
                                     parameter_estimate = NULL,
                                     prior_mean = NULL,
                                     zero_parameters = NULL,
                                     information_observed = NULL,
                                     information_expected = NULL,
                                     information_scale = c("total", "per_observation"),
                                     spbic_dimension = NULL) {
  scale_given <- !missing(information_scale)
  information_scale <- match.arg(information_scale)
  measures <- .sem_ic_match(measures, "measures")
  spbic_inputs <- c(parameter_estimate = !is.null(parameter_estimate),
                    prior_mean = !is.null(prior_mean),
                    zero_parameters = !is.null(zero_parameters),
                    information_observed = !is.null(information_observed),
                    spbic_dimension = !is.null(spbic_dimension))
  if (!"SPBIC" %in% measures && any(spbic_inputs)) {
    .dmar_refuse_argument(c(
      paste0("'", names(spbic_inputs)[spbic_inputs], "'", collapse = ", "),
      " ", if (sum(spbic_inputs) == 1L) "is" else "are", " used only for ",
      "SPBIC; add \"SPBIC\" to 'measures'."), arg = "measures")
  }
  if (!"IBIC" %in% measures && !is.null(information_expected)) {
    .dmar_refuse_argument(
      "'information_expected' is used only for IBIC; add \"IBIC\" to 'measures'.",
      arg = c("information_expected", "measures"))
  }
  pt <- NULL
  lik <- NULL
  identification <- NULL
  if (!is.null(fit)) {
    given <- c(log_likelihood = !is.null(log_likelihood), npar = !is.null(npar),
               N = !is.null(N),
               parameter_estimate = !is.null(parameter_estimate),
               information_observed = !is.null(information_observed),
               information_expected = !is.null(information_expected),
               information_scale = scale_given)
    if (any(given)) {
      .dmar_refuse_argument(c(
        "Supply either 'fit' or summary input, not both: a fitted model ",
        "supplies its own log likelihood, parameter count, sample size, and ",
        "per observation information matrices, which are multiplied by N (remove ",
        paste0("'", names(given)[given], "'", collapse = ", "), ")."),
        arg = "fit")
    }
    options <- .sem_fit_check(fit)
    lik <- .sem_likelihood_check(fit, options, "'fit'", measures)
    if (any(c("IBIC", "SPBIC") %in% measures)) {
      pt <- .sem_information_fit_check(fit, options, measures)
      identification <- .sem_identification(fit)
    }
    ll <- lavaan::logLik(fit)
    log_likelihood <- as.numeric(ll)
    npar <- as.numeric(attr(ll, "df"))
    N <- lavaan::lavInspect(fit, "ntotal")
    # A two-stage fit whose first stage read cases its N leaves out (cases
    # observed on an auxiliary variable alone): lavaan forms the log
    # likelihood as N times the complete data log likelihood per case at
    # the EM estimates, so its log likelihood and N rest on fewer cases
    # than the moments. Both are taken instead from the fit computed again
    # with every case the first stage read (.ts_all_rows() in
    # R/two_stage_internals.R). In one group, or in several that no
    # constraint holds together, the estimates stay, and each group's log
    # likelihood becomes its per case value times the cases read, with N
    # their total. In several groups held together by a constraint the
    # estimates move with the groups' counts, and the log likelihood is the
    # one at the estimates with every case read. Through the first versions
    # of 1.1.0 the criteria were lavaan's, with a warning; the warning
    # remains for a fit that cannot be computed again.
    ts_rows <- .ts_all_rows(fit)
    if (ts_rows$extra > 0L) {
      if (isTRUE(ts_rows$corrected)) {
        log_likelihood <- as.numeric(lavaan::logLik(ts_rows$fit))
        N <- sum(ts_rows$read)
      } else {
        .ts_aux_only_warning(
          fit,
          rest = paste0("the log likelihood, which lavaan forms as N times ",
                        "the complete data log likelihood per case, the ",
                        "sample size in the penalties of BIC, SABIC, and ",
                        "HBIC, and the N row rest"),
          consequence = paste0("every criterion, AIC included, weighs fit ",
                               "against complexity at a sample size other ",
                               "than the one the moments describe (DMAR ",
                               "could not compute the fit again with every ",
                               "case)"))
      }
    }
    if ("IBIC" %in% measures) {
      information_expected <- lavaan::lavInspect(fit, "information.expected")
    }
    if ("SPBIC" %in% measures) {
      parameter_estimate <- lavaan::coef(fit)
      information_observed <- lavaan::lavInspect(fit, "information.observed")
    }
    information_scale <- "per_observation"
    source <- paste0("lavaan ", utils::packageVersion("lavaan"))
  } else {
    if (is.null(log_likelihood) || is.null(npar) || is.null(N)) {
      .dmar_refuse_argument(
        "Supply 'fit', or the summary input 'log_likelihood', 'npar', and 'N'.",
        arg = c("fit", "log_likelihood", "npar"))
    }
    source <- "summary input"
  }
  .sem_scalar(log_likelihood, "log_likelihood")
  .sem_scalar(npar, "npar", 1, integer = TRUE)
  .sem_scalar(N, "N", 1, strict = TRUE)
  log_likelihood <- .sem_number(log_likelihood)
  N <- .sem_number(N)
  d <- as.integer(npar)
  if (!is.null(spbic_dimension)) {
    .sem_scalar(spbic_dimension, "spbic_dimension", 1, integer = TRUE)
    if (spbic_dimension > d) {
      .dmar_refuse_argument(c(
        "'spbic_dimension' must not exceed npar (", d, ")."),
        arg = "spbic_dimension")
    }
  }
  values <- stats::setNames(rep(NA_real_, length(.SEM_IC_NAMES)), .SEM_IC_NAMES)
  values["AIC"] <- -2 * log_likelihood + 2 * d
  values["BIC"] <- -2 * log_likelihood + d * log(N)
  values["SABIC"] <- -2 * log_likelihood + d * log((N + 2) / 24)
  values["HBIC"] <- -2 * log_likelihood + d * log(N / (2 * pi))
  if ("SABIC" %in% measures && N <= 22) {
    .dmar_warn(c(
      "With N = ", format(N), " the SABIC penalty per parameter, ",
      "log((N + 2)/24), is not positive, so SABIC does not penalize ",
      "complexity; read BIC or HBIC instead."), "result_qualified")
  }
  details <- list(source = source, N = N, npar = d)
  if (!is.null(lik)) {
    details <- c(details, list(likelihood = lik$likelihood, missing = lik$missing,
                               conditioning = lik$conditioning,
                               meanstructure = lik$meanstructure,
                               groups = lik$groups, estimator = lik$estimator))
    if (!is.null(identification)) details$identification <- identification
    if (isTRUE(ts_rows$extra > 0L)) {
      details$two_stage_rows <- ts_rows[c("read", "nobs", "extra",
                                          "corrected", "reestimated")]
    }
  }
  details$input_information_scale <- if (is.null(fit)) information_scale else
    "per observation (lavaan), multiplied by N"
  details$criterion_scale <- "-2 log likelihood plus a penalty; smaller is preferred"
  if ("IBIC" %in% measures) {
    if (is.null(information_expected)) {
      .dmar_refuse_argument(
        "IBIC needs 'information_expected', the expected information matrix.",
        arg = "information_expected")
    }
    ie <- .sem_information_matrix(information_expected, d, "information_expected",
                                  information_scale, N)
    values["IBIC"] <- -2 * log_likelihood - d * log(2 * pi) + ie$logdet
    details$expected_information <- if (is.null(fit)) "supplied" else
      paste0("lavaan's expected information at the fitted moments ",
             "(h1.information = \"structured\")")
    details$logdet_total_expected_information <- ie$logdet
  }
  if ("SPBIC" %in% measures) {
    if (!is.numeric(parameter_estimate) || length(parameter_estimate) != d ||
        any(!is.finite(parameter_estimate))) {
      .dmar_refuse_argument(c(
        "'parameter_estimate' must have one finite estimate per independent ",
        "free parameter (", d, ")."), arg = "parameter_estimate")
    }
    if (is.null(information_observed)) {
      .dmar_refuse_argument(
        "SPBIC needs 'information_observed', the observed information matrix.",
        arg = "information_observed")
    }
    io <- .sem_information_matrix(information_observed, d, "information_observed",
                                  information_scale, N,
                                  names(parameter_estimate))
    prior <- .sem_prior(parameter_estimate, prior_mean, zero_parameters, pt)
    # chol(I) = R with I = R'R, so Q = |R (estimate - prior)|^2, stable and
    # nonnegative by construction.
    Q <- sum((io$chol %*% (parameter_estimate - prior$mean))^2)
    if (!is.finite(Q)) .dmar_refuse_argument(
      "The SPBIC quadratic form overflowed.", arg = NULL)
    sd <- if (is.null(spbic_dimension)) d else as.integer(spbic_dimension)
    boundary <- Q <= sd
    penalty <- if (boundary) Q else sd * (1 + log(Q) - log(sd))
    if (boundary) {
      .dmar_warn(c(
        "SPBIC uses its boundary branch (Q = ", format(Q, digits = 4),
        " is at most the dimension ", sd, "), where the prior is a point ",
        "mass at the prior mean; its behavior there is little studied, so ",
        "interpret it cautiously (Bollen et al., 2012, Eq. 15; Bollen et ",
        "al., 2014, note 7)."), "result_qualified")
    }
    values["SPBIC"] <- -2 * log_likelihood + penalty
    details <- c(details, list(
      observed_information = if (is.null(fit)) "supplied" else
        "lavaan's negative Hessian of the log likelihood (observed.information = \"hessian\")",
      parameter_estimate = parameter_estimate, prior_mean = prior$mean,
      prior_rule = prior$rule,
      quadratic_form = Q, spbic_dimension = sd,
      dimension_rule = if (is.null(spbic_dimension)) "all independent free parameters"
                       else "supplied by spbic_dimension",
      penalty = penalty, boundary = boundary))
  }
  not_requested <- setdiff(.SEM_IC_NAMES, measures)
  values[not_requested] <- NA_real_
  if (any(!is.finite(values[measures]))) {
    .dmar_refuse_argument("A requested criterion is not finite.", arg = NULL)
  }
  details$footer <- .sem_information_footer(details, measures, not_requested)
  out <- .as_dmar_tbl(
    data.frame(term = c(.SEM_IC_NAMES, "log_likelihood", "npar", "N"),
               value = c(unname(values), log_likelihood, d, N),
               stringsAsFactors = FALSE),
    fixed_terms = c(.SEM_IC_NAMES, "log_likelihood"),
    subclass = "dmar_sem_information")
  attr(out, "sem_information") <- details
  out
}

# The lines written under a sem_information_criteria() table: the scale, the
# criteria not requested, the likelihood of a fitted model, and what IBIC
# and SPBIC depend on.
.sem_information_footer <- function(details, measures, not_requested) {
  out <- paste0("Smaller is preferred; only a difference between models fitted ",
                "to the same data under the same likelihood is meaningful.")
  if (length(not_requested)) {
    out <- c(out, paste0("Not requested: ", paste(not_requested, collapse = ", "),
                         " (see 'measures')."))
  }
  if (!is.null(details$likelihood)) {
    line <- paste0("Log likelihood: ", details$likelihood)
    if (length(details$conditioning)) {
      line <- paste0(line, ", conditional on ",
                     paste(details$conditioning, collapse = ", "))
    }
    out <- c(out, paste0(line, "."))
  }
  rows <- details$two_stage_rows
  if (isTRUE(rows$corrected)) {
    out <- c(out, paste0(
      "Two-stage fit: the log likelihood and N count the ", sum(rows$read),
      " cases the first stage read, ", rows$extra, " of them observed on an ",
      "auxiliary variable alone, which lavaan's logLik() and N (",
      sum(rows$nobs), ") leave out",
      if (isTRUE(rows$reestimated)) {
        paste0("; a constraint holds the groups together, so the ",
               "estimates were computed again with every case read, and the ",
               "log likelihood is at them")
      }, "."))
  }
  info <- intersect(c("IBIC", "SPBIC"), measures)
  if (length(info)) {
    id <- unique(details$identification)
    out <- c(out, paste0(
      paste(info, collapse = " and "), if (length(info) == 1L) " depends" else " depend",
      " on how the latent variables are identified",
      if (length(id)) paste0(" (here ", paste(id, collapse = " and "), ")") else "",
      " and on the mean structure", if ("IBIC" %in% info) ", IBIC also on the units" else "",
      "; compare only models that share them."))
  }
  if ("SPBIC" %in% measures) {
    out <- c(out, sprintf("SPBIC prior: %s; Q = %s, dimension %d%s.",
                          details$prior_rule, format(details$quadratic_form, digits = 6),
                          details$spbic_dimension,
                          if (isTRUE(details$boundary)) ", boundary branch" else ""))
  }
  out
}

# Every display of a sem_information_criteria() table (the console,
# knit_print(), as_kable()) writes these lines with it; see
# .dmar_tbl_footer() in R/dmar_output_helpers.R.
#' @exportS3Method .dmar_tbl_footer dmar_sem_information
#' @noRd
.dmar_tbl_footer.dmar_sem_information <- function(x, ...) {
  .dmar_footer_add(NextMethod(), attr(x, "sem_information")$footer,
                   "sem_information")
}
