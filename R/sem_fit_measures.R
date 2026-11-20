#' Fit Measures for a Structural Equation Model
#'
#' Reports in one table the measures researchers use to judge the global
#' fit of a structural equation model: the chi square test of exact fit,
#' the root mean square error of approximation (RMSEA) with its confidence
#' interval and the tests of close and not close fit, the incremental
#' indices CFI, TLI, NFI, IFI, and RNI, the SRMR, and gamma hat. The
#' function takes either a fitted \pkg{lavaan} model, whose own statistics,
#' sample size convention, and corrected (robust) versions it reads, or the
#' summary statistics a paper reports, so that a published fit table can
#' be reproduced and extended with an interval for the RMSEA and the tests
#' it implies. A fitted model requires \pkg{lavaan} (in Suggests).
#'
#' @param fit A fitted \pkg{lavaan} model. It must have converged to an
#'   admissible solution, as Details describes (lavaan's post-estimation
#'   check, made on the model's own parameters when the fit has auxiliary
#'   variables). Supply either \code{fit} or the summary statistics below,
#'   not both.
#' @param chisq,df For summary statistic input, the model's chi square test
#'   statistic \eqn{T} and its degrees of freedom: the likelihood ratio
#'   statistic of a normal theory maximum likelihood (ML) fit. Summary input
#'   gives the standard measures only.
#' @param N For summary statistic input, the total sample size, summed over
#'   groups: a whole number greater than \code{n_groups}.
#' @param baseline_chisq,baseline_df For summary statistic input, the chi
#'   square statistic and degrees of freedom of the baseline (independence)
#'   model fitted to the same data. The incremental indices need them and
#'   are \code{NA} without them. The baseline model is nested in the target
#'   model, so \code{baseline_df} must be positive and at least \code{df}.
#' @param p For summary statistic input, the number of observed variables
#'   (not the number of parameters). Gamma hat and adjusted gamma hat need
#'   it and are \code{NA} without it; they are defined here for a single
#'   group and a covariance structure (see Details), so they are \code{NA},
#'   with a warning, when \code{df} exceeds the \eqn{p(p + 1)/2} covariance
#'   moments, and \code{df} can be at most \eqn{p(p + 3)/2}. Below that
#'   bound summary input cannot show that a model restricts the means, and
#'   supplying \code{p} states that it does not.
#' @param n_groups For summary statistic input, the number of independent
#'   groups of a multiple group model, whose \code{chisq}, \code{df}, and
#'   \code{N} are totals over the groups. Default 1. A fitted model supplies
#'   its own number of groups.
#' @param conf_level Confidence level of the RMSEA interval and of the
#'   gamma hat and adjusted gamma hat intervals, which come from the same
#'   noncentrality limits. Defaults to 0.90, the conventional level of the
#'   RMSEA interval, whose two limits decide the tests of close and not
#'   close fit at .05 (see Details); \code{conf_level = 0.95} gives 95
#'   percent intervals.
#' @param close_fit The RMSEA value \eqn{\epsilon_0} under the null
#'   hypothesis of close fit, \eqn{H_0: \epsilon \le \epsilon_0}.
#'   Default 0.05.
#' @param not_close_fit The RMSEA value \eqn{\epsilon_1} under the null
#'   hypothesis of not close fit, \eqn{H_0: \epsilon \ge
#'   \epsilon_1}. Default 0.05.
#' @param version Which version of the measures to report for a fitted
#'   model: \code{"auto"} (the default), \code{"standard"},
#'   \code{"scaled"}, or \code{"robust"}. See Details. Summary statistic
#'   input accepts \code{"auto"} and \code{"standard"}, which coincide.
#' @param sample_size_convention For summary statistic input, the
#'   multiplier \eqn{a} in \eqn{T = a\hat{F}} with which \code{chisq} was
#'   formed: \code{"N-1"} (the default), each group's sample size less one,
#'   so \eqn{a = N - G} for \eqn{G} groups, or \code{"N"}, so \eqn{a = N}.
#'   lavaan's default normal likelihood and Mplus form \eqn{T} with
#'   \eqn{N}, so their statistics take \code{"N"}. A fitted model supplies
#'   its own convention.
#'
#' @details
#' \subsection{Exact Fit and the RMSEA}{
#' Write \eqn{T} for the model's chi square statistic on \eqn{df} degrees
#' of freedom, \eqn{a} for the sample size multiplier with which it was
#' formed (\eqn{T = a\hat{F}}, where \eqn{\hat{F}} is the minimized fit
#' function), \eqn{G} for the number of groups, and \eqn{\alpha} for 1
#' minus \code{conf_level}. The \emph{p}-value of the test of exact fit
#' (row \code{p_chi_square}) is the upper tail of the central chi square
#' distribution at \eqn{T}. Under a model that is wrong by an amount small
#' relative to sampling error, \eqn{T} is approximately noncentral chi
#' square with noncentrality \eqn{\lambda = a \cdot df \cdot \epsilon^2 /
#' G}, where \eqn{\epsilon} is the population RMSEA, the misfit per
#' degree of freedom on the scale of the fit function (Steiger & Lind,
#' 1980; Browne & Cudeck, 1993). The point estimate is \eqn{\hat{\epsilon}
#' = \sqrt{G \max(T - df, 0) / (a \cdot df)}}, and the confidence interval
#' inverts the noncentral chi square distribution at \eqn{T}
#' (\code{\link{ci_nc_chisq}} does the inversion) and maps its noncentrality
#' limits to the RMSEA metric in the same way. The factor \eqn{G} is
#' Steiger's (1998, Eqs. 23 and 24) multiple group RMSEA, which divides the
#' discrepancy by the average degrees of freedom per group; applying the
#' single group formula to a multiple group statistic gives a value too
#' small by a factor of \eqn{\sqrt{G}} (Steiger, 1998, Eq. 25).
#'
#' When \eqn{T} falls at or below the upper \eqn{\alpha/2} quantile of the
#' central chi square distribution, the lower limit is 0. When it falls at
#' or below the lower \eqn{\alpha/2} quantile, no noncentrality is
#' consistent with \eqn{T} at that level, and the upper limit is reported
#' as 0, as lavaan reports it and as Yuan, Chan, Marcoulides, and Bentler
#' (2016, Eq. 4) define the corresponding T-size RMSEA, so that the upper
#' limit and the test of not close fit always agree. This happens with
#' probability \eqn{\alpha/2} for a model that fits exactly, and more often
#' when the normal theory statistic is deflated, for example by data with
#' light tails.
#'
#' The RMSEA is conventionally reported with a 90 percent interval, the
#' one Steiger and Lind (1980) formed from the .05 and .95 percentile
#' points of the noncentrality, whose details Browne and Cudeck (1993)
#' give and which lavaan's summary prints. Its two limits are the two
#' one-sided tests at level .05 described below: a \eqn{100(1 - 2\alpha)}
#' percent interval settles both tests at level \eqn{\alpha}, so the lower
#' limit exceeds \eqn{\epsilon_0} exactly when close fit is rejected, and
#' the upper limit falls below \eqn{\epsilon_1} exactly when not close fit
#' is rejected (MacCallum, Browne, & Sugawara, 1996, p. 136 and Table 1,
#' with one bound playing both roles). A lower limit of zero likewise says
#' that the test of exact fit does not reject at .05. The default
#' \code{conf_level = 0.90} therefore gives the conventional interval,
#' which reproduces lavaan's printed limits, and every RMSEA interval in
#' DMAR defaults to it. \code{conf_level = 0.95}, the level of the
#' package's other intervals, gives a wider interval whose limits decide the
#' two tests at .025 rather than at .05. The printed table states the level
#' beneath it.
#'
#' A statistic in the millions, as from a large model fitted to a very
#' large sample, puts the noncentrality where R's noncentral chi square
#' distribution function, which lavaan uses, loses its accuracy (above
#' about \eqn{10^5}) and then fails (above about \eqn{2 \times 10^6}, where
#' it once returned an interval collapsed to a point below its estimate).
#' The limits and the tests of close and not close fit are then computed
#' here, from the \code{chi_square} and \code{df} rows, with the Poisson
#' mixture described in \code{\link{ci_nc_chisq}}; for a robust version
#' whose corrected statistic DMAR cannot rebuild (FIML, ordered data, or
#' two-stage estimation) they are \code{NA}, with a warning, once the
#' statistic exceeds \eqn{10^6}.
#' }
#'
#' \subsection{Tests of Close and Not Close Fit}{
#' The test of close fit (Browne & Cudeck, 1993) asks whether the data can
#' reject \eqn{H_0: \epsilon \le \epsilon_0}; its \emph{p}-value
#' is \eqn{P\{\chi^2_{df}(\lambda_0) \ge T\}} with \eqn{\lambda_0 = a \cdot
#' df \cdot \epsilon_0^2 / G} (row \code{p_close_fit}). Failing to
#' reject close fit does not establish it. The test of not close fit
#' (MacCallum, Browne, & Sugawara, 1996) reverses the roles: it asks
#' whether the data can reject \eqn{H_0: \epsilon \ge
#' \epsilon_1}, so that a small \emph{p}-value,
#' \eqn{P\{\chi^2_{df}(\lambda_1) \le T\}} with \eqn{\lambda_1 = a \cdot df
#' \cdot \epsilon_1^2 / G} (row \code{p_not_close_fit}), is evidence
#' that the population RMSEA is below \eqn{\epsilon_1}. It is the
#' equivalence test of \code{\link{equivalence_rmsea}}, which reports it
#' with its interval and verdict. MacCallum et al. (1996) used 0.05 for
#' both bounds, the defaults here; lavaan's default bound for not close fit
#' is 0.08, which \code{not_close_fit = 0.08} reproduces. These values are
#' conventions for framing a test, not cutoffs that certify a model.
#' }
#'
#' \subsection{Incremental Indices}{
#' The incremental indices compare the model with a baseline model, by
#' default the independence model, whose statistic and degrees of freedom
#' are \eqn{T_b} and \eqn{df_b}. The comparative fit index is \eqn{CFI = 1
#' - \max(T - df, 0) / \max(T - df, T_b - df_b, 0)}, equal to 1 when both
#' are zero (Bentler, 1990). Without the truncation at zero it is the
#' relative noncentrality index, \eqn{RNI = 1 - (T - df)/(T_b - df_b)}
#' (Bentler, 1990; McDonald & Marsh, 1990), which can fall outside
#' \eqn{[0, 1]}. The Tucker-Lewis index (Tucker & Lewis, 1973), also called
#' the nonnormed fit index (Bentler & Bonett, 1980), is \eqn{TLI =
#' (T_b/df_b - T/df)/(T_b/df_b - 1)}; it rewards parsimony and can also
#' exceed 1. The normed fit index is \eqn{NFI = (T_b - T)/T_b} (Bentler &
#' Bonett, 1980), and the incremental fit index is \eqn{IFI = (T_b -
#' T)/(T_b - df)} (Bollen, 1989, Eq. 3). No index is clipped to \eqn{[0,
#' 1]}, and an index whose denominator is zero or negative is \code{NA}.
#' Brosseau-Liard and Savalei (2014, Eqs. 3 to 10) give the two population
#' values these sample indices estimate: the NFI, CFI, RNI, and IFI all
#' converge to the population CFI, \eqn{1 - F_H/F_B} for the population
#' ML discrepancies of the model and the baseline (the NFI poorly, since
#' it does not subtract the degrees of freedom), and the TLI to the
#' population TLI, \eqn{1 - (F_H/F_B)(df_B/df_H)}. None of them has a
#' confidence interval here:
#' the usual asymptotic theory assumes small misfit, which the independence
#' model violates.
#'
#' A fitted model with no degrees of freedom reproduces the moments
#' exactly, so its statistic is zero in every version, and so is the model
#' term \eqn{T - c \cdot df} of the robust indices, whatever the scaling
#' constant \eqn{c}, which is 0/0 there (lavaan's robust and scaled indices
#' of such a model come out \code{NA} or not by rounding). Its CFI is
#' therefore 1, and its RNI, NFI, and IFI are their values at a statistic
#' of zero, those of lavaan's standard version, in every version that has
#' them; its TLI, like the RMSEA, divides by the degrees of freedom and is
#' \code{NA}. A model with no degrees of freedom need not fit exactly: an
#' inequality constraint that binds, such as an error variance held at
#' zero by \code{bounds = "pos.var"}, fixes its parameter as an equality
#' would, though lavaan does not count it among the degrees of freedom.
#' Such a model's statistic is not zero (32.8 for one factor of three
#' variables correlated .8, .8, and .5, \eqn{N} = 200), and its CFI, RNI,
#' NFI, and IFI are formed from it (.925, .925, .926, and .926), with the
#' standard NFI computed as \eqn{1 - T/T_b} where lavaan reports 1 at zero
#' degrees of freedom; the TLI and the RMSEA rows are \code{NA} as before,
#' and the warning and the footer say that the model does not fit exactly.
#' }
#'
#' \subsection{SRMR and Gamma Hat}{
#' The standardized root mean square residual (row \code{srmr}) is the
#' root mean square of the covariance residuals, each divided by the
#' product of the two sample standard deviations (lavaan's
#' \code{srmr_bentler_nomean}). It is read from a fitted model and is
#' \code{NA} for summary input. It is the same quantity in every version,
#' since no test statistic enters it. lavaan's own \code{srmr} also
#' averages in the mean residuals when the model has a mean structure,
#' which makes it smaller by about \eqn{\sqrt{t/(t + p)}}, for \eqn{t = p(p
#' + 1)/2} covariance moments of \eqn{p} variables, when the means are
#' saturated, as under FIML; this table keeps the covariance definition
#' throughout. With auxiliary variables (\pkg{lavaan}'s \code{aux} under
#' FIML) the SRMR is over the model's own variables, the auxiliary
#' variables set aside as in the admissibility check below. \pkg{lavaan}'s
#' measures average over them as well, and an auxiliary variable, a
#' saturated correlate, adds residuals of about zero that lower the SRMR of
#' the same model: by about 8 percent with one auxiliary variable and 15
#' percent with two for a three-factor model of nine items. The test
#' statistic, and the measures built from it, compare the model with the
#' saturated model with the auxiliary variables in both, so they need no
#' such change. For the unbiased estimate of the SRMR, its interval, and its
#' tests see \code{\link{sem_residual_fit}}.
#'
#' Gamma hat rescales the estimated noncentrality to a goodness of fit
#' metric, \eqn{\hat{\gamma} = p/[p + 2(T - df)/a]} for \eqn{p} observed
#' variables, and adjusted gamma hat penalizes it for the degrees of
#' freedom used, \eqn{1 - [p(p + 1)/(2 df)](1 - \hat{\gamma})}, counting
#' the \eqn{p(p + 1)/2} covariance moments of one group. They estimate the
#' population values of the goodness of fit index (GFI) and the adjusted
#' GFI that LISREL reports, which Steiger (1989) showed to be these
#' functions of the population noncentrality for a model invariant under a
#' constant scaling factor (a result he restates in Steiger, 1990, p.
#' 178). Neither is truncated, so both exceed 1 when \eqn{T < df}. Each is
#' a monotone function of the noncentrality, so the noncentrality limits of
#' the RMSEA interval give their limits, \eqn{p/(p + 2\lambda/a)} at each
#' limit \eqn{\lambda} and the adjusted index of that (Steiger & Fouladi,
#' 1997, pp. 251--252), in the rows \code{gamma_hat_ci_lower} through
#' \code{adjusted_gamma_hat_ci_upper}. The limits are at most 1, the upper
#' one is 1 when the lower RMSEA limit is 0, and a point estimate above 1
#' lies above its interval.
#'
#' This gamma hat is the index lavaan 0.7-3 reports as \code{gfi_lrt},
#' with the same limits; lavaan caps the point estimate at 1 and reports it
#' only when it is requested by name. The \code{gfi} lavaan prints by
#' default is the same function of Browne's residual based (reweighted
#' least squares) statistic in place of \eqn{T}, which differs from gamma
#' hat in a finite sample (0.959 against 0.957 for the three-factor model
#' of the examples), and lavaan's \code{gfi_lisrel} and \code{agfi_lisrel}
#' are the sample GFI and adjusted GFI of LISREL, not estimates of their
#' population values. lavaan also reports robust versions, and so does this
#' table: for a robust ML fit of complete continuous data with a
#' Satorra-Bentler, Yuan-Bentler, scaled and shifted, or mean and variance
#' adjusted statistic (MLM, MLR, MLMV, MLMVS), the robust version
#' estimates the noncentrality as \eqn{s(T_s - df_s)}, for the corrected
#' statistic \eqn{T_s}, its scaling constant \eqn{s}, and its degrees of
#' freedom \eqn{df_s}, which is \eqn{T - c \cdot df} in each case (see the
#' robust version below), the estimate behind the robust RMSEA of
#' Brosseau-Liard, Savalei, and Li (2012), and takes its limits from the
#' robust RMSEA interval; for the first three this is lavaan's
#' \code{gfi_lrt.robust} and its limits.
#'
#' Steiger's result is for covariance structures. Neither index is computed
#' for a fitted model that restricts the means, whose statistic and degrees
#' of freedom include the misfit and the restrictions of the means while
#' \eqn{p} and \eqn{p(p + 1)/2} count the covariances: a growth model, or
#' intercepts constrained equal (by a shared label, with or without
#' lavaan's \code{ceq.simple = TRUE}), constrained to a constant, entered
#' in a linear constraint, tied to a loading, or held at the bound of an
#' inequality constraint (which lavaan's degrees of freedom do not count).
#' Whether a fit restricts its means is read from the derivatives of its
#' model implied moments, with the constraints that bind at the estimates
#' imposed (the equality constraints, and the inequality constraints
#' lavaan treats as active when it computes the standard errors): the
#' means are free when the parameter directions that leave every
#' covariance unchanged move the means in as many independent directions
#' as there are modeled means, so that the means add no degree of freedom
#' (the means of exogenous
#' variables under \code{fixed.x = TRUE} are fixed at their sample values
#' and are not counted). A mean structure whose means are free, as
#' lavaan's default parameterization leaves them or as one indicator's
#' intercept fixed at zero with the latent mean free leaves them, changes
#' neither index. Summary input carries no mean structure to read, so it
#' shows restricted means only when \code{df} exceeds \eqn{p(p + 1)/2},
#' which only restricted means can produce, and the indices are then
#' \code{NA}; below that bound they are computed from the summary
#' statistics on the understanding that the model is a covariance
#' structure (the statistic of a linear growth model of four waves, 12.2
#' on 5 df with \eqn{N} = 400, gives 0.991 and 0.982, which the fitted
#' model does not). The two indices are computed for the standard and
#' robust versions just described of a single group, single level ML fit
#' of continuous variables without FIML or sampling weights, and for
#' summary input when \code{p} is given, and are \code{NA} otherwise.
#' }
#'
#' \subsection{Standard, Scaled, and Robust Versions}{
#' A fitted model can carry more than one test statistic, and lavaan
#' computes the fit measures from each.
#' \describe{
#' \item{\code{"standard"}}{The measures above computed from the
#'   uncorrected statistic. They are the right ones for a normal theory ML
#'   fit of complete continuous data.}
#' \item{\code{"scaled"}}{The same formulas with a scaled (or mean and
#'   variance adjusted) statistic \eqn{T_s} on its degrees of freedom
#'   \eqn{df_s} put in place of \eqn{T} and \eqn{df}: the RMSEA and its
#'   interval invert the noncentral chi square at \eqn{T_s} on
#'   \eqn{df_s}, the \code{chi_square} and \code{df} rows of the same
#'   table (Brosseau-Liard, Savalei, & Li, 2012, Eqs. 11 and 13, where
#'   \eqn{df_s = df} for a scaled statistic; Savalei, 2018, Table 1), and
#'   the tests of close and not close fit, which those sources do not
#'   treat, follow the same substitution so that they keep their duality
#'   with the interval. This is the substitution EQS and Mplus report, and
#'   the one lavaan reported when Brosseau-Liard et al. (2012, p. 911)
#'   wrote, and the one a published table of scaled measures was computed
#'   with. It is not
#'   lavaan's default for a Satorra-Bentler or Yuan-Bentler statistic (MLM,
#'   MLR) or a mean and variance adjusted one (MLMVS): lavaan 0.7-3 refers
#'   the unscaled \eqn{T} to the noncentral chi square on \eqn{c \cdot df}
#'   degrees of freedom, the trace of \eqn{U\Gamma} for scaling constant
#'   \eqn{c}, which gives the same point estimate but a different interval
#'   and different close fit \emph{p}-values from those the rows beside them
#'   imply. A fit made with lavaan's \code{rmsea.scaled.mplus = TRUE} gives
#'   the values reported here, as lavaan's default does for a scaled and
#'   shifted statistic (MLMV, WLSMV). The scaled RMSEA and CFI do not
#'   consistently estimate the population values the standard formulas
#'   define (Brosseau-Liard et al., 2012; Brosseau-Liard & Savalei, 2014;
#'   Savalei, 2018, 2021); they are offered to reproduce a published
#'   table.}
#' \item{\code{"robust"}}{lavaan's corrected measures, which estimate the
#'   population values of the standard measures: for a robust ML fit of
#'   complete continuous data (MLM, MLR), the nonnormality corrections of
#'   Brosseau-Liard, Savalei, and Li (2012) for the RMSEA and of
#'   Brosseau-Liard and Savalei (2014) for the CFI, TLI, and RNI, with the
#'   correction of Savalei (2018) for a scaled and shifted statistic
#'   (MLMV). The point estimates do not depend on which corrected
#'   statistic accompanies them: each subtracts from the ML statistic
#'   \eqn{T} its asymptotic expectation \eqn{c \cdot df}, with \eqn{c} the
#'   trace of \eqn{U\Gamma} over \eqn{df} (Brosseau-Liard et al., 2012,
#'   Appendix; Savalei, 2018, Eq. 6). lavaan computes no robust measures
#'   beside the older mean and variance adjusted statistic (MLMVS,
#'   \code{test = "mean.var.adjusted"}, the statistic Savalei, 2018,
#'   Appendix B, writes \eqn{T_{MVold}} and no longer recommends; her mean
#'   and variance adjusted statistic is the scaled and shifted one), whose
#'   scaling constant \eqn{s} and degrees of freedom \eqn{df_s} give
#'   \eqn{c = s \cdot df_s / df}, so they are computed here, for the model and for its
#'   baseline, with the interval and the tests of close and not close fit
#'   from \eqn{T/c} on \eqn{df} (Brosseau-Liard et al., 2012, Eq. 10); they
#'   equal lavaan's robust measures for the same model fitted by MLM, and
#'   the footer says that DMAR computed them (for a saturated model, whose
#'   \eqn{c} is 0/0 and whose robust CFI and RNI are 1 at a statistic of
#'   zero, nothing is computed from \eqn{c}, and the footer names the
#'   corrections alone, as for MLM). For ordered categorical data
#'   (WLSMV), the correction of
#'   Savalei (2021); and for a FIML fit (\code{missing = "ml"}), the FIML
#'   corrected indices of Zhang and Savalei (2023). The last are not
#'   corrections for nonnormality: they evaluate the complete data ML fit
#'   function at the FIML estimates, with the saturated FIML estimates in
#'   place of the sample moments, and so approximate the RMSEA and CFI the
#'   model would have had with complete data, whereas the uncorrected FIML
#'   indices change with the amount of missing data. lavaan uses their
#'   version 3, which does not assume normality, when the fit has a scaled
#'   statistic (MLR), and their version 6, which does, otherwise. The
#'   corrected values rest on the corrected statistic, which under FIML
#'   differs from the chi square shown in the table. After two-stage
#'   estimation they are the corrected indices of Zhang and Savalei (2023),
#'   which correct the complete data indices of the second stage for the
#'   missing data: after robust two-stage estimation (\code{missing =
#'   "robust.two.stage"}) their TS V2, the version they recommend, and after
#'   normal theory two-stage estimation (\code{"two.stage"}) its normal
#'   theory analog, which takes the asymptotic covariance matrix of the
#'   stage one moments to be the inverse of their observed information where
#'   TS V2 takes the sandwich of the observed and first order information,
#'   and which they did not study (see \code{\link{fa_confirmatory}}). Under FIML and
#'   two-stage estimation alike, Zhang and Savalei define only the RMSEA and
#'   the CFI; the robust TLI and RNI, the RMSEA interval, and the tests of
#'   close and not close fit are lavaan's extension of the correction. NFI
#'   and IFI have no robust version, in lavaan or here, and are \code{NA};
#'   gamma hat has the one described above.}
#' }
#' The default \code{"auto"} reports the robust version for a fit that
#' carries a corrected test statistic (MLM, MLR, WLSMV, a sampling weighted
#' fit with lavaan's default test, and the like) or uses FIML, and the
#' standard version otherwise. A fit to ordered data by DWLS or ULS has only
#' the standard statistic unless a corrected test is requested (it is the
#' default test of \code{estimator = "DWLS"} and \code{"ULS"}, which
#' \code{"WLSMV"} and \code{"ULSMV"} extend with a corrected one), and
#' lavaan gives that statistic no reference distribution, so \code{"auto"}
#' reports its standard measures with a warning that the statistic is not
#' chi square distributed and that a corrected test is needed for the
#' robust measures. The same holds for a DWLS fit of continuous data whose
#' test was set to \code{"standard"} (lavaan's default test for continuous
#' data adds Browne's residual based statistic, which is chi square under
#' normality and is then the one reported; the warning says to leave out
#' \code{test = "standard"} or to use a corrected test), for a ULS fit of
#' continuous data whose test was set so (see \emph{ULS Fits of Continuous
#' Data}), and for a sampling weighted fit whose test was set to
#' \code{"standard"}, whose statistic is a weighted pseudo likelihood ratio
#' statistic. Full WLS is the exception among the
#' estimators for ordered data: its weight matrix is the inverse of the
#' asymptotic covariance matrix of the thresholds and polychoric
#' correlations, its statistic is asymptotically chi square (Browne, 1984;
#' Muthén, 1984), lavaan reports its \emph{p}-value, and \code{"auto"}
#' reports its standard measures without a warning (a request for the
#' robust version says that they are the ones to report). The reason in
#' the footer names the estimator. The measures of a ULS fit of continuous
#' data are reported only where they do not depend on the units of the
#' variables (see \emph{ULS Fits of Continuous Data}). Browne's residual
#' based statistic, which lavaan computes beside the standard one, is an
#' alternative test and not
#' a correction, and does not count. When lavaan provides a corrected
#' statistic but no robust measures for it (and DMAR cannot compute them),
#' \code{"auto"} reports, with a warning, the standard version for a
#' continuous ML fit, whose standard measures estimate the population
#' values of the ML measures as the robust ones do and the scaled ones do
#' not (Savalei, 2018), though their interval and tests assume normal
#' data, and the scaled version for any other fit. The \code{chi_square}, \code{df}, and
#' \code{p_chi_square} rows show the scaled statistic under the scaled and
#' robust versions when the fit has one, and the standard statistic
#' otherwise. The \code{sem_fit} attribute and the printed footer name the
#' version reported, why, and the statistic used.
#' }
#'
#' \subsection{ULS Fits of Continuous Data}{
#' The corrections of Brosseau-Liard and colleagues target the population
#' discrepancy of the ML fit function. Their use with another fit
#' function needs separate justification (Savalei, 2018). For continuous
#' ULS, \code{version = "robust"} is therefore refused for the entire
#' robust family, including CFI, TLI, and RNI. Invariance to multiplying
#' every variable by the same constant is not enough to establish the
#' population interpretation of an incremental index. This is a
#' conservative supported-method policy, not a claim that Savalei proved
#' every continuous-ULS index invalid. Ordered-data corrections with their
#' own derivation are handled separately.
#'
#' With ULSM or ULSMV, \code{version = "auto"} selects eligible scaled
#' indices. These put a corrected test statistic into the index formulas;
#' they are descriptive approximations, do not estimate ML-population
#' fit, and have no established license to use conventional ML cutoffs.
#' When the fit restricts means, incompatible units in the mean and
#' covariance portions of the ULS discrepancy can make these indices
#' change even under a common rescaling. Those rows are withheld, with
#' the existing mean-structure qualification.
#'
#' The standard statistic \eqn{N F_{ULS}} is in squared-covariance units.
#' Its RMSEA and incremental indices that subtract degrees of freedom
#' change with those units and are \code{NA}; the standard NFI, a ratio
#' of two such statistics, remains available when eligible. This does
#' not apply to a separately justified residual-based statistic. A fit
#' by \code{estimator = "ULS"} alone is tested in
#' lavaan 0.7-3 by Browne's normal theory residual based statistic, which
#' weights the residuals by the inverse of their normal theory covariance
#' matrix and has no units as long as the estimates at which it is
#' evaluated do not change with the units either. They do not when the
#' model leaves the means free, or restricts them through parameters the
#' covariances do not share, as the equal intercepts above do, and the
#' standard measures formed from that statistic are then reported as they
#' are. When the model restricts the means through parameters the
#' covariances share, as the loadings enter the means through the latent
#' means of a multiple group model whose loadings and intercepts are equal
#' (the scalar rung of \code{\link{measurement_invariance}}), the ULS fit
#' function trades the mean residuals against the covariance residuals,
#' which are in other units, so the estimates change with the units, and
#' so do Browne's statistic and every measure formed from it: on the scalar
#' rung of the first example of \code{\link{measurement_invariance}}, the
#' CFI is .8874 for the raw scores and .8861 for the scores divided by 10,
#' and the RMSEA .0680 and .0684. The RMSEA, its interval, the tests of
#' close and not close fit, and the CFI, TLI, RNI, NFI, and IFI are then
#' \code{NA}, with the same warning, element, and line; the
#' \code{chi_square} row, whose chi square reference distribution holds at
#' any root-\eqn{n} consistent estimates (Browne, 1984, Proposition 4), and
#' the SRMR, which describe the fit at hand, are reported. Whether the two parts
#' share a parameter is read from the same derivatives as whether the
#' means are restricted: they share one when the rank of the derivatives of
#' every model implied moment falls short of the rank of the covariance
#' rows plus that of the mean rows. With ordered data the robust
#' measures of a fit with \code{"ULSMV"} are the categorical correction
#' of Savalei (2021), formed from polychoric correlations, which have no
#' units, and none of this applies. Maximum likelihood is scale invariant;
#' \code{estimator = "MLR"} is the estimator for nonnormal continuous data.
#' }
#'
#' \subsection{Two-Stage Fits}{
#' A fit by two-stage maximum likelihood (Yuan & Bentler, 2000; Savalei &
#' Bentler, 2009) fits the model to the EM estimates of the means and
#' covariances as if they were the moments of complete data. lavaan's
#' standard statistic for it is that naive one, which is not chi square
#' distributed, so its standard measures (TS V0 in the terms of Zhang &
#' Savalei, 2023) are refused, and so are its scaled measures, which put an
#' approximate statistic into the standard formulas; the default
#' \code{"auto"} reports the robust version. The \code{chi_square} row is
#' then the residual based statistic of Savalei and Bentler (2009, Eq. 8; Eq. 14 with auxiliaries)
#' after normal theory two-stage estimation (\code{missing = "two.stage"}),
#' computed by DMAR from lavaan's moments and the observed information of
#' the saturated model, and lavaan's scaled statistic after robust two-stage
#' estimation (Savalei & Falk, 2014): the statistics
#' \code{\link{fa_confirmatory}} and \code{\link{measurement_invariance}} report for
#' the same fits. The indices are not computed from either statistic but
#' from the naive one and its expectation, as Zhang and Savalei (2023)
#' define them: TS V2 after robust two-stage estimation and its normal
#' theory analog after normal theory two-stage estimation.
#'
#' A case observed on an auxiliary variable (lavaan's \code{aux}) but on
#' none of the model's variables is read by the first stage, whose EM
#' estimates of the means and covariances it moves, and is left out of the
#' sample size lavaan 0.7-3 reports (\code{lavaan::lavInspect(fit,
#' "nobs")}), which counts it as empty. The naive statistic, and with it
#' lavaan's scaled statistic, the corrected indices, the RMSEA multiplier,
#' and the \code{N} row, then rest on fewer cases than the moments, while
#' the residual based statistic of a normal theory two-stage fit uses every
#' case, so the rows of the table do not describe one set of cases (and
#' lavaan's two-stage standard errors are too large by the square root of
#' the ratio of the two counts). \code{sem_fit_measures()} reports the
#' rows of a fit as \pkg{lavaan} computed them, and warns when a fit has
#' such cases and counts them. \code{\link{fa_confirmatory}(missing = "two_stage")}
#' and \code{\link{measurement_invariance}} keep such cases and compute
#' their fits again with the sample size set to every case the first stage
#' read, so that the fits they return, and the tables of this function for
#' them, count every case; FIML (\code{missing = "ml"}) with the same
#' auxiliary variables keeps them as well. Setting the cases aside instead
#' would condition the analysis on an observed variable of the model,
#' which biases the estimates when the items went missing together at
#' random given the auxiliary variables.
#' }
#'
#' \subsection{The Sample Size Convention}{
#' The multiplier \eqn{a} must be the one with which \eqn{T} was formed. A
#' fitted model supplies it: lavaan forms \eqn{T = N\hat{F}} for ML under
#' its default normal likelihood, and \eqn{T = (N - G)\hat{F}} under
#' \code{likelihood = "wishart"} and for the least squares estimators, and
#' the table follows the fit, so its values equal lavaan's (for the scaled
#' version, those lavaan gives under \code{rmsea.scaled.mplus = TRUE}; see
#' above). For summary
#' input the default \code{"N-1"} is the convention of Browne and Cudeck
#' (1993), MacCallum et al. (1996), Kelley and Lai (2011), and DMAR's
#' \code{\link{ci_rmsea}}, \code{\link{ss_aipe_rmsea}}, and
#' \code{\link{ss_power_sem}}. A statistic reported by lavaan under its
#' default likelihood or by Mplus was formed with \eqn{N} and takes
#' \code{sample_size_convention = "N"}. The wrong multiplier changes the
#' RMSEA, its limits, the close fit tests, and gamma hat by an amount of
#' order \eqn{1/N}.
#' }
#'
#' \subsection{What Is Refused}{
#' A fitted model is refused when it did not converge (an error of class
#' \code{dmar_not_converged}) and when its solution is inadmissible (a
#' refusal of class \code{dmar_inadmissible}): fit measures of such a
#' solution describe estimates outside the parameter space and do not
#' speak to the specified model. A solution is inadmissible when
#' \pkg{lavaan}'s post-estimation check finds a negative variance estimate
#' (a Heywood case among them) or a covariance matrix that is not positive
#' definite: that of the latent variables, that of the residuals of the
#' observed variables (theta), or that of the residuals of the structural
#' part of the model (psi). The message names the variable or the matrix,
#' and the condition carries them (see \code{\link{dmar_metadata}}). With
#' auxiliary variables (lavaan's \code{aux} under FIML) the check is made
#' on the model's own parameters, with the auxiliary variables set aside:
#' lavaan models an auxiliary variable as a saturated correlate, whose
#' covariances with the items run through their residuals, and its own
#' check of theta, which includes those covariances, fails once an
#' auxiliary variable is related to a factor strongly enough, however
#' proper the model's own solution. For one factor whose items have
#' uncorrelated residuals, and an auxiliary variable related to the items
#' only through the factor, that happens when the squared correlation
#' between the auxiliary variable and the factor exceeds
#' \eqn{(1 - H) / H}, with \emph{H} the maximal reliability of the
#' factor's items (coefficient \emph{H}; see \code{\link{reliability_H}}),
#' which is a correlation above .33 when \emph{H} is .9. Auxiliary
#' variables are recognized when they are given through lavaan's
#' \code{aux} argument; saturated correlates written into the model by
#' hand, as \code{semTools::cfa.auxiliary()} writes them, stay in theta
#' and are judged with it, so refit with \code{aux} to have them set
#' aside. The same check serves \code{\link{sem_residual_fit}},
#' \code{\link{sem_information_criteria}}, \code{\link{compare_sem}},
#' \code{\link{ecvi}}, and \code{\link{sem_exact_fit_bootstrap}}, and
#' \code{\link{ci_bifactor_indices}} judges its fit and its replications
#' by it.
#' Resolve the improper solution first, by changing the model or by
#' judging that a negative variance is within sampling error of zero and
#' fixing it, or pass the reported statistics as summary input if you
#' accept that responsibility. Multilevel fits are refused, since one set
#' of indices mixes the two levels; lavaan's level specific fit measures
#' (\code{lavaan::fitMeasures(fit, level = ...)}) are the tool for them. A
#' fit whose cases were removed by listwise deletion is evaluated on the
#' complete cases with a warning, since that assumes the data are missing
#' completely at random; refit it with \code{missing = "ml"} to use every
#' case.
#'
#' For a fit by normal theory two-stage estimation (\code{missing =
#' "two.stage"}) the residual based statistic, and with it the table, is
#' refused when the fit has sampling weights or clustered cases, conditions
#' on an observed exogenous variable (\code{fixed.x = TRUE}, the default of
#' \code{lavaan::sem()} for a model that regresses on one, or
#' \code{conditional.x = TRUE}), has inequality constraints or a parameter
#' on a bound, or has auxiliary variables that the installed lavaan does
#' not keep where DMAR reads them. Each message says what to do instead;
#' robust two-stage estimation, whose scaled statistic lavaan computes, is
#' the alternative throughout. A fit that conditions on an exogenous
#' variable is refused with an error of class
#' \code{dmar_two_stage_fixed_x}, which carries the exogenous variables,
#' and the advice to refit with \code{fixed.x = FALSE}, and with
#' \code{conditional.x = FALSE} too, since \pkg{lavaan} does not take the
#' first beside \code{conditional.x = TRUE}. The refit leaves the
#' estimates as they are when the exogenous variables are complete, and
#' otherwise also keeps the cases \pkg{lavaan} deletes for a missing
#' exogenous variable when it conditions on them (robust two-stage
#' estimation deletes them too), so that its estimates differ; to keep the
#' fit's estimates, refit on the rows the fit used, those on which every
#' exogenous variable is observed, which
#' \code{lavaan::lavInspect(fit, "case.idx")} lists.
#' }
#'
#' \subsection{Reading the Measures}{
#' No row of the table certifies a model, and the table reports no
#' cutoff. The conventions in use are descriptions of a population RMSEA:
#' MacCallum et al. (1996, p. 134) credit the view that values below .05
#' indicate close fit to Steiger (1989) and others, describe .05 to .08
#' as fair and above .10 as poor after Browne and Cudeck (1993), and call
#' .08 to .10 mediocre themselves. Chen, Curran, Bollen, Kirby, and
#' Paxton (2008) found no empirical support for .05, or any other value,
#' as a universal cutoff: a fixed bound rejects too many correct models
#' in samples of 100 or fewer, too few misspecified ones in samples of
#' 800 or more, and the test of close fit with a fixed bound performs
#' inconsistently across sample sizes, degrees of freedom, and
#' specifications. They recommend reading the RMSEA beside the chi square
#' statistic and the other indices rather than as a single accept or reject
#' decision; the interval and the two tests above are how this table
#' states the uncertainty, and the researcher chooses the bound.
#'
#' The cutoffs in common use for the other indices, near .95 for the CFI
#' and TLI, .06 for the RMSEA, and .08 for the SRMR, applied in pairs,
#' come from the simulations of Hu and Bentler (1999), who themselves
#' wrote that no single cutoff works equally well across conditions and
#' that the trade-off between the two kinds of error at samples of 250 or
#' fewer cannot be judged from the sample. Marsh, Hau, and Wen (2004)
#' showed that the rationale behind those simulations treated a minor
#' misspecification as a model to be rejected, so that the rules reward
#' rejecting models that fit acceptably and cannot hold across sample
#' sizes; they warn against reading the rules of thumb as golden rules.
#' What a value of an index means depends on the model: on the size of
#' the loadings and unique variances, which lower the chi square, RMSEA,
#' and SRMR of a misspecified model when the loadings are low (Heene,
#' Hilbert, Draxler, Ziegler, & Buehner, 2011); on the number of observed
#' variables, which inflates the statistic and every index formed from it
#' (Moshagen, 2012); and on the reliability of the factors, the number of
#' items and factors, and the degrees of freedom, which is why McNeish and
#' Wolf (2023) derive cutoffs for the model at hand by simulation, at
#' several sizes of misspecification, so that an index reads as an effect
#' size rather than as a verdict. This table reports the indices and their
#' uncertainty; the bound a model is held to is a judgment about the
#' model, stated in advance.
#'
#' With few degrees of freedom the RMSEA is unstable (Kenny, Kaniskan, &
#' McCoach, 2015): for correctly specified models with 1 to 3 degrees of
#' freedom and \eqn{N} = 400, the sample RMSEA exceeded .05 in 13 to 15
#' percent of samples, against none at 50 degrees of freedom, and its
#' interval at \eqn{N} = 50 and 1 degree of freedom was more than
#' sixteen times as wide as at \eqn{N} = 1,000 and 100 degrees of
#' freedom. They advise against computing the RMSEA when the degrees of
#' freedom are small and suggest diagnosing poor fit by fitting the
#' models that free the omitted parameters instead
#' (\code{\link{compare_sem}} compares them). The interval in this table
#' makes the instability visible rather than hiding it.
#'
#' The incremental indices measure fit against the baseline model, and
#' their metric is the baseline's misfit: under maximum likelihood the
#' independence baseline's statistic is \eqn{-(N - 1)\log|R|} for the
#' sample correlation matrix \eqn{R} (van Laar & Braeken, 2022). Two
#' correctly specified models with the same chi square and RMSEA
#' distributions have different CFI distributions when their
#' \eqn{|R|} differ, and misspecified models can rank differently on the
#' CFI and on the chi square, so a CFI is comparable across data sets only
#' when the baselines are; van Laar and Braeken recommend reporting
#' \eqn{|R|} and the baseline statistic and reading an incremental index
#' as the share of the baseline's misfit the model removes. A fitted
#' model's baseline statistic and degrees of freedom are
#' \code{lavaan::fitMeasures(fit, c("baseline.chisq", "baseline.df"))};
#' summary input takes them as given.
#' }
#'
#' The raw-data H1 screen is sufficient and nonexhaustive. In models
#' with fixed or conditional observed predictors, it screens response
#' variables only, since predictor moments are not free. It does not
#' assess conditional design rank or exact response-on-predictor relations.
#' No detected obstacle is not a proof of likelihood existence or uniqueness.
#'
#' @return A \code{dmar_tbl} (class \code{dmar_sem_fit}) with columns
#'   \code{term} and \code{value} and, always in this order, the rows
#'   \code{chi_square}, \code{df}, \code{p_chi_square} (the test of exact
#'   fit); \code{rmsea}, \code{rmsea_ci_lower}, \code{rmsea_ci_upper};
#'   \code{close_fit_bound}, \code{p_close_fit}, \code{not_close_fit_bound},
#'   \code{p_not_close_fit}; \code{cfi}, \code{tli}, \code{nfi},
#'   \code{ifi}, \code{rni}; \code{srmr}; \code{gamma_hat},
#'   \code{gamma_hat_ci_lower}, \code{gamma_hat_ci_upper},
#'   \code{adjusted_gamma_hat}, \code{adjusted_gamma_hat_ci_lower},
#'   \code{adjusted_gamma_hat_ci_upper}; \code{N} and \code{n_groups}. A row the input
#'   or the version does not compute is \code{NA}. The \code{sem_fit}
#'   attribute records the source, the version requested and reported and
#'   why, the statistic in the \code{chi_square} row, the correction behind
#'   the robust measures, the sample size convention and multiplier, why
#'   rows of a ULS fit of continuous data are \code{NA} (\code{uls};
#'   \code{NA} otherwise), the boundary conventions that applied, and the
#'   lines of the printed footer.
#'
#' @references
#' Bentler, P. M. (1990). Comparative fit indexes in structural models.
#'   \emph{Psychological Bulletin, 107}(2), 238--246.
#'   \doi{10.1037/0033-2909.107.2.238}
#'
#' Bentler, P. M., & Bonett, D. G. (1980). Significance tests and goodness
#'   of fit in the analysis of covariance structures.
#'   \emph{Psychological Bulletin, 88}, 588--606.
#'
#' Bollen, K. A. (1989). A new incremental fit index for general structural
#'   equation models. \emph{Sociological Methods & Research, 17}(3),
#'   303--316. \doi{10.1177/0049124189017003004}
#'
#' Brosseau-Liard, P. E., & Savalei, V. (2014). Adjusting incremental fit
#'   indices for nonnormality. \emph{Multivariate Behavioral Research,
#'   49}(5), 460--470. \doi{10.1080/00273171.2014.933697}
#'
#' Brosseau-Liard, P. E., Savalei, V., & Li, L. (2012). An investigation of
#'   the sample performance of two nonnormality corrections for RMSEA.
#'   \emph{Multivariate Behavioral Research, 47}(6), 904--930.
#'   \doi{10.1080/00273171.2012.715252}
#'
#' Browne, M. W. (1984). Asymptotically distribution-free methods for the
#'   analysis of covariance structures. \emph{British Journal of
#'   Mathematical and Statistical Psychology, 37}(1), 62--83.
#'
#' Browne, M. W., & Cudeck, R. (1993). Alternative ways of assessing
#'   model fit. In K. A. Bollen & J. S. Long (Eds.), \emph{Testing
#'   structural equation models} (pp. 136--162). Sage.
#'
#' Chen, F., Curran, P. J., Bollen, K. A., Kirby, J., & Paxton, P. (2008).
#'   An empirical evaluation of the use of fixed cutoff points in RMSEA
#'   test statistic in structural equation models. \emph{Sociological
#'   Methods & Research, 36}(4), 462--494. \doi{10.1177/0049124108314720}
#'
#' Heene, M., Hilbert, S., Draxler, C., Ziegler, M., & Buehner, M. (2011).
#'   Masking misfit in confirmatory factor analysis by increasing unique
#'   variances: A cautionary note on the usefulness of cutoff values of
#'   fit indices. \emph{Psychological Methods, 16}(3), 319--336.
#'   \doi{10.1037/a0024917}
#'
#' Hu, L., & Bentler, P. M. (1999). Cutoff criteria for fit indexes in
#'   covariance structure analysis: Conventional criteria versus new
#'   alternatives. \emph{Structural Equation Modeling, 6}(1), 1--55.
#'   \doi{10.1080/10705519909540118}
#'
#' Kelley, K., & Lai, K. (2011). Accuracy in parameter estimation for the
#'   root mean square error of approximation: Sample size planning for
#'   narrow confidence intervals.
#'   \emph{Multivariate Behavioral Research, 46}(1), 1--32.
#'   \doi{10.1080/00273171.2011.543027}
#'
#' Kenny, D. A., Kaniskan, B., & McCoach, D. B. (2015). The performance of
#'   RMSEA in models with small degrees of freedom. \emph{Sociological
#'   Methods & Research, 44}(3), 486--507. \doi{10.1177/0049124114543236}
#'
#' MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
#'   analysis and determination of sample size for covariance structure
#'   modeling. \emph{Psychological Methods, 1}(2), 130--149.
#'   \doi{10.1037/1082-989X.1.2.130}
#'
#' Marsh, H. W., Hau, K.-T., & Wen, Z. (2004). In search of golden rules:
#'   Comment on hypothesis-testing approaches to setting cutoff values for
#'   fit indexes and dangers in overgeneralizing Hu and Bentler's (1999)
#'   findings. \emph{Structural Equation Modeling, 11}(3), 320--341.
#'   \doi{10.1207/s15328007sem1103_2}
#'
#' McDonald, R. P., & Marsh, H. W. (1990). Choosing a multivariate model:
#'   Noncentrality and goodness of fit. \emph{Psychological Bulletin, 107},
#'   247--255.
#'
#' McNeish, D., & Wolf, M. G. (2023). Dynamic fit index cutoffs for
#'   confirmatory factor analysis models. \emph{Psychological Methods,
#'   28}(1), 61--88. \doi{10.1037/met0000425}
#'
#' Moshagen, M. (2012). The model size effect in SEM: Inflated
#'   goodness-of-fit statistics are due to the size of the covariance
#'   matrix. \emph{Structural Equation Modeling, 19}(1), 86--98.
#'   \doi{10.1080/10705511.2012.634724}
#'
#' Muthén, B. (1984). A general structural equation model with dichotomous,
#'   ordered categorical, and continuous latent variable indicators.
#'   \emph{Psychometrika, 49}(1), 115--132.
#'
#' Savalei, V. (2018). On the computation of the RMSEA and CFI from the
#'   mean-and-variance corrected test statistic with nonnormal data in SEM.
#'   \emph{Multivariate Behavioral Research, 53}(3), 419--429.
#'   \doi{10.1080/00273171.2018.1455142}
#'
#' Savalei, V. (2021). Improving fit indices in structural equation
#'   modeling with categorical data. \emph{Multivariate Behavioral
#'   Research, 56}(3), 390--407. \doi{10.1080/00273171.2020.1717922}
#'
#' Savalei, V., & Bentler, P. M. (2009). A two-stage approach to missing
#'   data: Theory and application to auxiliary variables. \emph{Structural
#'   Equation Modeling, 16}, 477--497. \doi{10.1080/10705510903008238}
#'
#' Savalei, V., & Falk, C. F. (2014). Robust two-stage approach outperforms
#'   robust full information maximum likelihood with incomplete nonnormal
#'   data. \emph{Structural Equation Modeling, 21}(2), 280--302.
#'   \doi{10.1080/10705511.2014.882692}
#'
#' Steiger, J. H. (1989). \emph{EzPATH: A supplementary module for SYSTAT
#'   and SYGRAPH}. SYSTAT.
#'
#' Steiger, J. H. (1990). Structural model evaluation and modification: An
#'   interval estimation approach. \emph{Multivariate Behavioral Research,
#'   25}(2), 173--180.
#'
#' Steiger, J. H. (1998). A note on multiple sample extensions of the RMSEA
#'   fit index. \emph{Structural Equation Modeling, 5}(4), 411--419.
#'   \doi{10.1080/10705519809540115}
#'
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
#'   estimation and the evaluation of statistical models. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' Steiger, J. H., & Lind, J. C. (1980). \emph{Statistically-based tests
#'   for the number of common factors}. Paper presented at the annual Spring
#'   meeting of the Psychometric Society, Iowa City, IA.
#'
#' Tucker, L. R., & Lewis, C. (1973). A reliability coefficient for maximum
#'   likelihood factor analysis. \emph{Psychometrika, 38}(1), 1--10.
#'
#' van Laar, S., & Braeken, J. (2022). Caught off base: A note on the
#'   interpretation of incremental fit indices. \emph{Structural Equation
#'   Modeling, 29}(6), 935--943. \doi{10.1080/10705511.2022.2050730}
#'
#' Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods for
#'   mean and covariance structure analysis with nonnormal missing data.
#'   \emph{Sociological Methodology, 30}, 165--200.
#'
#' Yuan, K.-H., Chan, W., Marcoulides, G. A., & Bentler, P. M. (2016).
#'   Assessing structural equation models by equivalence testing with
#'   adjusted fit indexes. \emph{Structural Equation Modeling, 23}(3),
#'   319--330. \doi{10.1080/10705511.2015.1065414}
#'
#' Zhang, X., & Savalei, V. (2023). New computations for RMSEA and CFI
#'   following FIML and TS estimation with missing data.
#'   \emph{Psychological Methods, 28}(2), 263--283.
#'   \doi{10.1037/met0000445}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{sem_residual_fit}} for the unbiased SRMR with its
#'   interval and tests, \code{\link{equivalence_rmsea}} for the test of not
#'   close fit with its verdict, \code{\link{ci_rmsea}} for an interval from
#'   a reported RMSEA, \code{\link{sem_information_criteria}} and
#'   \code{\link{compare_sem}} for comparing models, and the vignette
#'   \code{vignette("sem_fit_measures", package = "DMAR")}.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' # The statistics lavaan reports for the three-factor model of nine of
#' # the holzinger_swineford tests (N = 301), with its baseline.
#' # lavaan forms the statistic with N under its default likelihood.
#' sem_fit_measures(chisq = 85.3055, df = 24, N = 301,
#'                  baseline_chisq = 918.8516, baseline_df = 36, p = 9,
#'                  sample_size_convention = "N")
#'
#' # The same model fitted here. The interval is the conventional 90
#' # percent one, the default, whose limits are the two one-sided .05 tests
#' # of close and not close fit.
#' fit <- lavaan::cfa(
#'   "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
#'    textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
#'    speed   =~ t10_addition + t12_counting_groups_of_dots +
#'               t13_straight_and_curved_capitals",
#'   data = holzinger_swineford, std.lv = TRUE)
#' sem_fit_measures(fit)
#'
#' # The 95 percent interval, wider, decides the two tests at .025.
#' sem_fit_measures(fit, conf_level = 0.95)
#'
#' # A robust ML fit reports the robust measures by default.
#' fit_mlr <- lavaan::cfa(
#'   "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
#'    textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
#'    speed   =~ t10_addition + t12_counting_groups_of_dots +
#'               t13_straight_and_curved_capitals",
#'   data = holzinger_swineford, std.lv = TRUE, estimator = "MLR")
#' sem_fit_measures(fit_mlr)
#'
#' # The model fitted in each of the two schools: the multiple group RMSEA
#' # needs the number of groups (Steiger, 1998).
#' sem_fit_measures(chisq = 115.8513, df = 48, N = 301, n_groups = 2,
#'                  sample_size_convention = "N")
#'
#' @export
sem_fit_measures <- function(fit = NULL, chisq = NULL, df = NULL, N = NULL,
                             baseline_chisq = NULL, baseline_df = NULL,
                             p = NULL, n_groups = 1, conf_level = 0.90,
                             close_fit = 0.05, not_close_fit = 0.05,
                             version = c("auto", "standard", "scaled", "robust"),
                             sample_size_convention = c("N-1", "N")) {
  version <- match.arg(version)
  .sem_probability(conf_level, "conf_level")
  .sem_scalar(close_fit, "close_fit", 0, strict = TRUE)
  .sem_scalar(not_close_fit, "not_close_fit", 0, strict = TRUE)
  close_fit <- .sem_number(close_fit)
  not_close_fit <- .sem_number(not_close_fit)
  if (!is.null(fit)) {
    given <- c(chisq = !is.null(chisq), df = !is.null(df), N = !is.null(N),
               baseline_chisq = !is.null(baseline_chisq),
               baseline_df = !is.null(baseline_df), p = !is.null(p),
               n_groups = !missing(n_groups),
               sample_size_convention = !missing(sample_size_convention))
    if (any(given)) {
      .dmar_refuse_argument(c(
        "Supply either 'fit' or the summary statistics, not both: a fitted ",
        "model supplies its own statistics, number of groups, and sample size ",
        "convention (remove ",
        paste0("'", names(given)[given], "'", collapse = ", "), ")."),
        arg = "fit")
    }
    res <- .sem_fit_measures_lavaan(fit, conf_level, close_fit, not_close_fit,
                                    version)
  } else {
    sample_size_convention <- match.arg(sample_size_convention)
    res <- .sem_fit_measures_summary(chisq, df, N, baseline_chisq, baseline_df,
                                     p, n_groups, conf_level, close_fit,
                                     not_close_fit, version,
                                     sample_size_convention)
  }
  out <- .as_dmar_tbl(data.frame(term = names(res$values),
                                 value = unname(res$values),
                                 stringsAsFactors = FALSE),
                      conf_level = conf_level,
                      p_terms = c("p_chi_square", "p_close_fit", "p_not_close_fit"),
                      subclass = "dmar_sem_fit")
  res$details$footer <- c(res$details$footer, res$details$boundary)
  attr(out, "sem_fit") <- res$details
  out
}

# The rows of every sem_fit_measures() table, in order.
.SEM_FIT_TERMS <- c("chi_square", "df", "p_chi_square",
                    "rmsea", "rmsea_ci_lower", "rmsea_ci_upper",
                    "close_fit_bound", "p_close_fit",
                    "not_close_fit_bound", "p_not_close_fit",
                    "cfi", "tli", "nfi", "ifi", "rni", "srmr",
                    "gamma_hat", "gamma_hat_ci_lower", "gamma_hat_ci_upper",
                    "adjusted_gamma_hat", "adjusted_gamma_hat_ci_lower",
                    "adjusted_gamma_hat_ci_upper", "N", "n_groups")

# Gamma hat and adjusted gamma hat for one group, with their limits (see the
# help page). `excess` is the estimated population discrepancy, (T - df)/a
# for the standard version and c (T_s - df)/a for the robust one, not
# truncated at zero; the limits are the same function of the noncentrality
# limits of the version's RMSEA interval, df RMSEA^2 = lambda/a for one
# group (Steiger & Fouladi, 1997, pp. 251--252), and adjusted gamma hat is
# an increasing function of gamma hat, so its limits follow.
.sem_gamma_hat <- function(excess, df, p, rmsea_lower, rmsea_upper) {
  gamma <- function(discrepancy) {
    den <- p + 2 * discrepancy
    if (is.finite(den) && den > 0) p / den else NA_real_
  }
  adjust <- function(g) {
    if (df > 0 && is.finite(g)) 1 - p * (p + 1) / (2 * df) * (1 - g) else NA_real_
  }
  g <- c(gamma(excess), gamma(df * rmsea_upper^2), gamma(df * rmsea_lower^2))
  c(gamma_hat = g[1L], gamma_hat_ci_lower = g[2L], gamma_hat_ci_upper = g[3L],
    adjusted_gamma_hat = adjust(g[1L]), adjusted_gamma_hat_ci_lower = adjust(g[2L]),
    adjusted_gamma_hat_ci_upper = adjust(g[3L]))
}

# A ratio that is NA when its denominator is not positive.
.sem_ratio <- function(num, den) if (is.finite(den) && den > 0) num / den else NA_real_

# Why "auto" reports the standard measures of a fit with no corrected
# statistic, said from lavaan's options rather than from a choice the
# user is assumed to have made. The standard test is lavaan's default for
# estimator = "DWLS", "ULS", and "WLS" with ordered data, as for the
# continuous ML, GLS, ULS, and WLS estimators; with sampling weights its
# default is a corrected test, so a weighted fit without one had its test
# set to "standard". Full WLS is the asymptotically distribution free
# estimator of Browne (1984) for continuous data and its categorical
# analog (Muthén, 1984), whose statistic needs no correction, and the
# least squares fits of continuous data are tested in lavaan 0.7-3 by
# Browne's normal theory residual based statistic (`standard_test`)
# unless their test was set to "standard", whose statistic lavaan gives
# no reference distribution, as it does for DWLS and ULS with ordered
# data.
.sem_standard_reason <- function(estimator, categorical, weighted,
                                 standard_test, correlation = FALSE) {
  if (weighted) {
    return(paste0("a sampling weighted fit whose test was set to ",
                  "\"standard\", so it has no corrected statistic"))
  }
  if (categorical) {
    if (identical(estimator, "WLS")) {
      return(paste0("a fit to ordered data by full WLS, whose statistic is ",
                    "asymptotically chi square without a correction"))
    }
    return(paste0("a fit to ordered data by ", estimator, " with the ",
                  "standard test alone, so it has no corrected statistic",
                  switch(estimator,
                         DWLS = " (estimator = \"WLSMV\" is DWLS with one)",
                         ULS = " (estimator = \"ULSMV\" is ULS with one)",
                         "")))
  }
  if (estimator %in% c("DWLS", "ULS", "PML") &&
      identical(standard_test, "standard")) {
    return(paste0("a ", estimator, " fit of complete data with the standard ",
                  "test alone, whose statistic lavaan gives no reference ",
                  "distribution",
                  if (identical(estimator, "ULS")) {
                    if (correlation) {
                      " and which is the ULS discrepancy of a correlation structure"
                    } else {
                      " and which is in the units of the squared covariances"
                    }
                  }))
  }
  switch(estimator,
    ML = "a normal theory ML fit of complete data with no corrected statistic",
    GLS = "a normal theory GLS fit of complete data with no corrected statistic",
    WLS = paste0("an asymptotically distribution free (WLS) fit of complete ",
                 "data, whose statistic needs no correction"),
    paste0("a ", estimator, " fit of complete data with no corrected statistic",
           if (startsWith(standard_test, "browne.residual")) {
             ", tested by Browne's residual based statistic"
           }))
}

# The summary statistic path: the closed forms from T, df, N, and the
# baseline statistics.
.sem_fit_measures_summary <- function(chisq, df, N, baseline_chisq, baseline_df,
                                      p, n_groups, conf_level, close_fit,
                                      not_close_fit, version, convention) {
  if (version %in% c("scaled", "robust")) {
    .dmar_refuse_argument(c(
      "Summary statistics give the standard measures only: the scaled and ",
      "robust versions need the scaling constants and corrected statistics of ",
      "a fitted model. Supply 'fit', or use version = \"standard\"."),
      arg = "fit")
  }
  if (is.null(chisq) || is.null(df) || is.null(N)) {
    .dmar_refuse_argument(
      "Supply 'fit', or the summary statistics 'chisq', 'df', and 'N'.",
      arg = c("fit", "chisq", "df"))
  }
  .sem_scalar(chisq, "chisq", 0)
  .sem_scalar(df, "df", 0, integer = TRUE)
  .sem_scalar(n_groups, "n_groups", 1, integer = TRUE)
  .sem_scalar(N, "N", 1, strict = TRUE, integer = TRUE)
  chisq <- .sem_number(chisq)
  df <- .sem_number(df)
  n_groups <- .sem_number(n_groups)
  N <- .sem_number(N)
  if (N <= n_groups) {
    .dmar_refuse_argument(
      "'N' is the total sample size over the groups and must exceed 'n_groups'.",
      arg = c("N", "n_groups"))
  }
  if (xor(is.null(baseline_chisq), is.null(baseline_df))) {
    .dmar_refuse_argument(
      "Supply both 'baseline_chisq' and 'baseline_df', or neither.",
      arg = c("baseline_chisq", "baseline_df"))
  }
  has_baseline <- !is.null(baseline_chisq)
  if (has_baseline) {
    .sem_scalar(baseline_chisq, "baseline_chisq", 0)
    .sem_scalar(baseline_df, "baseline_df", 0, strict = TRUE, integer = TRUE)
    baseline_chisq <- .sem_number(baseline_chisq)
    baseline_df <- .sem_number(baseline_df)
    if (baseline_df < df) {
      .dmar_refuse_argument(c(
        "The baseline model is nested in the target model, so 'baseline_df' ",
        "must be at least 'df'."), arg = c("baseline_df", "df"))
    }
  }
  if (!is.null(p)) {
    .sem_scalar(p, "p", 1, integer = TRUE)
    p <- .sem_number(p)
    if (n_groups > 1) {
      .dmar_refuse_argument(c(
        "'p' is used only for gamma hat, which is defined here for a single ",
        "group; leave 'p' out when 'n_groups' is greater than 1."),
        arg = c("p", "n_groups"))
    }
    if (df > p * (p + 3) / 2) {
      .dmar_refuse_argument(c(
        "'df' exceeds the number of sample moments of p observed variables, ",
        "p(p + 3)/2 = ", p * (p + 3) / 2, "."), arg = "df")
    }
  }
  a <- .sem_multiplier(N, n_groups, convention)
  vals <- stats::setNames(rep(NA_real_, length(.SEM_FIT_TERMS)), .SEM_FIT_TERMS)
  vals[c("chi_square", "df", "N", "n_groups")] <- c(chisq, df, N, n_groups)
  vals[c("close_fit_bound", "not_close_fit_bound")] <- c(close_fit, not_close_fit)
  boundary <- character(0)
  if (df > 0) {
    vals[["p_chi_square"]] <- stats::pchisq(chisq, df, lower.tail = FALSE)
    rows <- .sem_rmsea_rows(chisq, df, a, n_groups, conf_level, close_fit,
                            not_close_fit)
    vals[names(rows)] <- rows
    if (stats::pchisq(chisq, df) <= (1 - conf_level) / 2) {
      boundary <- .sem_boundary_text("no_upper_root")
    }
  }
  if (has_baseline) {
    num <- max(chisq - df, 0)
    den <- max(chisq - df, baseline_chisq - baseline_df, 0)
    vals[["cfi"]] <- if (den == 0) 1 else 1 - num / den
    if (df > 0) {
      vals[["tli"]] <- .sem_ratio(baseline_chisq / baseline_df - chisq / df,
                                  baseline_chisq / baseline_df - 1)
    }
    vals[["nfi"]] <- .sem_ratio(baseline_chisq - chisq, baseline_chisq)
    vals[["ifi"]] <- .sem_ratio(baseline_chisq - chisq, baseline_chisq - df)
    rni_den <- baseline_chisq - baseline_df
    vals[["rni"]] <- if (rni_den != 0) 1 - (chisq - df) / rni_den else NA_real_
    incremental <- c("tli", "nfi", "ifi", "rni")
    undefined <- incremental[is.na(vals[incremental])]
    if (df == 0) undefined <- setdiff(undefined, "tli")
    if (length(undefined)) {
      .dmar_warn(c(
        "These baseline statistics leave ",
        paste(undefined, collapse = ", "),
        " undefined (a zero or negative denominator), so ",
        if (length(undefined) > 1L) "they are" else "it is", " NA."),
        "result_qualified")
    }
  }
  gamma_note <- if (is.null(p)) "not computed: 'p' was not supplied" else "computed"
  if (!is.null(p)) {
    if (df > p * (p + 1) / 2) {
      gamma_note <- paste0("not computed: df exceeds the p(p + 1)/2 covariance ",
                           "moments, so the model restricts the means")
      .dmar_warn(c(
        "'df' (", format(df), ") exceeds the p(p + 1)/2 = ", p * (p + 1) / 2,
        " covariances of p = ", format(p), " variables, so the model ",
        "restricts the means; gamma hat and adjusted gamma hat are defined ",
        "for covariance structures and are NA."), "result_qualified")
    } else {
      gamma <- .sem_gamma_hat((chisq - df) / a, df, p, vals[["rmsea_ci_lower"]],
                              vals[["rmsea_ci_upper"]])
      vals[names(gamma)] <- gamma
    }
  }
  if (df == 0) {
    exact <- .sem_exact_fit(c(chisq = chisq, df = df,
                              baseline.chisq = if (has_baseline) baseline_chisq))
    if (exact) {
      vals[c("rmsea", "rmsea_ci_lower", "rmsea_ci_upper")] <- 0
      if (has_baseline) vals[["tli"]] <- 1
    }
    boundary <- .sem_boundary_text(if (exact) "zero_df" else "zero_df_misfit")
    .sem_warn_zero_df(exact)
  }
  details <- list(
    source = "summary statistics",
    version_requested = version, version = "standard",
    version_reason = "summary statistics give the standard measures only",
    test = "the supplied chi square statistic",
    correction = NA_character_,
    sample_size_convention = convention, multiplier = a,
    N = N, groups = n_groups,
    srmr = "not available from summary statistics",
    gamma = gamma_note,
    uls = NA_character_,
    boundary = boundary)
  details$footer <- c(
    "Version: standard (summary statistics give the standard measures only)",
    sprintf("Sample size: N = %s in %s; multiplier a = %s (sample_size_convention = \"%s\")",
            format(N), if (n_groups == 1) "1 group" else paste(n_groups, "groups"),
            format(a), convention))
  list(values = vals, details = details)
}

# The fitted model path: lavaan's own values under the version reported,
# with the fit's own sample size convention and number of groups.
.sem_fit_measures_lavaan <- function(fit, conf_level, close_fit, not_close_fit,
                                     version_requested) {
  options <- .sem_fit_check(fit)
  ts_rows <- .ts_all_rows(fit)
  if (isTRUE(ts_rows$extra > 0L)) {
    if (!isTRUE(ts_rows$corrected)) {
      .dmar_refuse(c("Two-stage fit measures require the same all-row sample ",
                     "as the first-stage moments; the auxiliary-only case ",
                     "correction could not be computed. Retain FIML or refit ",
                     "a supported two-stage model."),
                   "two_stage_rows_unavailable",
                   rows = ts_rows[c("read", "nobs", "extra")])
    }
    fit <- ts_rows$fit
    options <- .sem_fit_check(fit)
  }
  h1_diagnostic <- .fiml_fit_h1_diagnostic(fit)
  if (!is.null(h1_diagnostic)) {
    .dmar_refuse(c("Fit indices requiring the unrestricted H1 model are ",
                   "unavailable. ", h1_diagnostic, " The constrained model's ",
                   "parameter estimates and H0-based standard errors may still ",
                   "be usable; inspect those separately rather than reporting ",
                   "these fit indices."),
                 "h1_unavailable", diagnostic = h1_diagnostic)
  }

  if (lavaan::lavInspect(fit, "nlevels") > 1L) {
    .dmar_refuse_argument(c(
      "sem_fit_measures() does not evaluate multilevel fits: one set of ",
      "indices mixes the within and between levels. Use lavaan's level ",
      "specific fit measures, lavaan::fitMeasures(fit, level = ...)."),
      arg = NULL)
  }
  # R's pchisq() warns that it did not converge each time lavaan evaluates
  # it at a noncentrality above about 2e6 (dozens of times for one interval).
  # Every limit and close fit test lavaan forms there is recomputed or set
  # to NA below, so those warnings are muffled; any other warning passes.
  fm <- withCallingHandlers(
    lavaan::fitMeasures(fit, "all",
                        fm.args = list(rmsea.ci.level = conf_level,
                                       rmsea.h0.closefit = close_fit,
                                       rmsea.h0.notclosefit = not_close_fit)),
    warning = function(cond) {
      if (grepl("pnchisq", conditionMessage(cond), fixed = TRUE))
        invokeRestart("muffleWarning")
    })
  if (!all(c("chisq", "df") %in% names(fm)) || !is.finite(fm[["chisq"]])) {
    .dmar_refuse_argument(c(
      "The fit has no test statistic (it was fitted with test = \"none\"); ",
      "refit it with a test."), arg = NULL)
  }
  conv <- .sem_lavaan_convention(fit)
  categorical <- length(lavaan::lavNames(fit, "ov.ord")) > 0L
  corrected <- .lavaan_corrected_tests(options$test)
  fiml <- options$missing %in% c("ml", "ml.x")
  two_stage <- options$missing %in% c("two.stage", "robust.two.stage")
  weighted <- isTRUE(options$.sampling.weights)
  scaled_test <- attr(fm, "scaled.test")
  if (is.null(scaled_test) && length(corrected)) scaled_test <- corrected[[1L]]
  standard_test <- attr(fm, "standard.test")
  if (is.null(standard_test)) standard_test <- "standard"
  # The robust measures of a continuous ML fit with the older mean and
  # variance adjusted statistic, which lavaan does not compute, computed
  # here at the requested level and bounds (.sem_robust_mva()).
  mva <- .sem_robust_mva(fit, conf_level, close_fit, not_close_fit)
  if (length(mva)) fm <- c(fm, mva)
  # A continuous ML fit, complete or listwise: its robust measures estimate
  # the population values of the ML measures, and so do its standard ones,
  # with the small sample bias corrected as if the data were normal, but
  # its scaled ones do not (Savalei, 2018).
  continuous_ml <- identical(options$estimator, "ML") && !categorical &&
    !fiml && !two_stage
  has_scaled <- "chisq.scaled" %in% names(fm)
  # A ULS fit of continuous data: lavaan's robust RMSEA, and the standard
  # one when the standard statistic is (N - G) times the ULS fit function,
  # change with the units of the variables (.sem_uls_continuous()), so
  # those rows are NA below.
  uls_continuous <- .sem_uls_continuous(fit)
  uls_units <- uls_continuous && identical(standard_test, "standard")
  has_robust <- all(c("rmsea.robust", "cfi.robust") %in% names(fm))
  if (uls_continuous && version_requested == "robust") {
    .dmar_refuse_argument(c(
      "The ML-derived robust fit-index correction is not supported for ",
      "continuous ULS fits. Invariance of an incremental index to a common ",
      "change of units does not establish its population interpretation ",
      "under a different fit function (Savalei, 2018). Use version = ",
      "\"scaled\" with a corrected ULS test for qualified scaled indices, ",
      "or refit by ML/MLR for ML-population fit indices."), arg = NULL)
  }
  if (uls_continuous) has_robust <- FALSE
  # Whether the standard statistic lacks a chi square reference
  # distribution. lavaan 0.7-3's lav_model_test() gives the standard
  # statistic of DWLS, ULS, and PML no reference distribution and no
  # p-value, whatever the data, and that of full WLS, whose weight matrix
  # is the inverse of the asymptotic covariance matrix of the moments (of
  # the thresholds and polychoric correlations with ordered data), the
  # chi square distribution it has asymptotically (Browne, 1984; Muthén,
  # 1984), with its p-value. With continuous data the default test of
  # DWLS and ULS adds Browne's residual based statistic, which is chi
  # square under normality and which fitMeasures() then reports as the
  # standard one (`standard_test`), so only a test set to "standard"
  # leaves the statistic without a reference distribution.
  not_chisq <- isTRUE(options$estimator %in% c("DWLS", "ULS", "PML")) &&
    identical(standard_test, "standard")

  # A two-stage fit has two statistics lavaan reports beside its robust
  # measures, and neither is the test to report with them. Its standard
  # statistic treats the EM moments as the moments of complete data, so it
  # is not chi square distributed and its measures (Zhang and Savalei's TS
  # V0) are biased in small samples. Its scaled measures put an approximate
  # statistic into the standard formulas, which changes the population
  # values they estimate (Brosseau-Liard, Savalei, & Li, 2012;
  # Brosseau-Liard & Savalei, 2014). Both are refused in favor of the
  # corrected indices of Zhang and Savalei (2023): after robust two-stage
  # estimation their TS V2, which they recommend, and after normal theory
  # two-stage estimation its normal theory analog (see .ts_indices_label()).
  if (two_stage && version_requested %in% c("standard", "scaled")) {
    .dmar_refuse_argument(c(
      "version = \"", version_requested, "\" is not available for a ",
      "two-stage fit (missing = \"", options$missing, "\"): ",
      if (version_requested == "standard") {
        paste0("its standard statistic treats the EM estimates of the ",
               "means and covariances as complete data, so it is not chi ",
               "square distributed")
      } else {
        paste0("its scaled measures put an approximate statistic into ",
               "the standard formulas, which changes the population values ",
               "they estimate")
      },
      ". Use version = \"robust\" (the default \"auto\" chooses it), ",
      "the corrected indices of Zhang and Savalei (2023)."), arg = NULL)
  }

  # Which version to report, and why.
  if (version_requested == "auto" && uls_continuous && has_scaled) {
    version <- "scaled"
    reason <- paste0("continuous ULS: scaled descriptive indices from the ",
                     "corrected test; ML-derived robust indices are unsupported")
  } else if (version_requested == "auto") {
    if (length(corrected) || fiml) {
      why <- if (fiml) {
        "the fit uses FIML (missing = \"ml\")"
      } else if (two_stage) {
        paste0("the fit uses two-stage estimation (missing = \"",
               options$missing, "\")")
      } else {
        paste0("the fit carries the corrected statistic ", scaled_test)
      }
      if (has_robust) {
        version <- "robust"
        reason <- why
      } else if (has_scaled && continuous_ml) {
        version <- "standard"
        reason <- paste0(why, ", for which lavaan computes no robust measures")
        .dmar_warn(c(
          "lavaan computes no robust fit measures for the ", scaled_test,
          " statistic, so the standard versions are reported; they ",
          "estimate the population values of the ML measures, as the ",
          "robust ones do and the scaled ones do not (Savalei, 2018), but ",
          "their p-value, RMSEA interval, and close fit tests assume ",
          "normal data. Refit with estimator = \"MLM\", \"MLR\", or ",
          "\"MLMV\" for the robust measures."), "result_qualified")
      } else if (has_scaled) {
        version <- "scaled"
        reason <- paste0(why, ", for which lavaan computes no robust measures")
        .dmar_warn(c(
          "lavaan computes no robust fit measures for the ", scaled_test,
          " statistic, so the scaled versions are reported; they put the ",
          "scaled statistic in the standard formulas and do not estimate ",
          "the population RMSEA and CFI consistently."), "result_qualified")
      } else {
        version <- "standard"
        reason <- paste0(why, ", but lavaan reports no corrected measures for it")
      }
    } else {
      version <- "standard"
      reason <- .sem_standard_reason(options$estimator, categorical,
                                     weighted, standard_test,
                                     correlation = isTRUE(options$correlation))
    }
  } else {
    version <- version_requested
    reason <- "requested"
    if (version == "robust" && !has_robust) {
      .dmar_refuse_argument(c(
        "lavaan computes robust fit measures only for a fit with a corrected ",
        "test statistic (for example estimator = \"MLR\" or \"MLM\", or WLSMV ",
        "for ordered data) or a FIML fit (missing = \"ml\"); this fit has ",
        "neither, ",
        if (not_chisq && !categorical) {
          paste0("and its standard statistic, that of ", options$estimator,
                 ", is not chi square distributed (lavaan gives it no ",
                 "p-value). Leave out test = \"standard\", so that lavaan ",
                 "tests the fit by Browne's residual based statistic, or ",
                 "refit it with a corrected test (estimator = ",
                 if (identical(options$estimator, "ULS")) {
                   "\"ULSM\" or \"ULSMV\""
                 } else "\"WLSMV\" with ordered = FALSE", ").")
        } else if (not_chisq || weighted) {
          paste0("and its standard statistic is not chi square distributed. ",
                 "Refit it with a corrected test (",
                 if (categorical) "estimator = \"WLSMV\", lavaan's default for ordered data"
                 else "estimator = \"MLR\", or leave out test = \"standard\" so that lavaan's default for sampling weighted data applies",
                 ").")
        } else if (categorical) {
          paste0("and its standard statistic, that of full WLS, is ",
                 "asymptotically chi square (Browne, 1984; Muth\u00e9n, 1984), ",
                 "so its standard measures are the ones to report. The ",
                 "corrected measures of Savalei (2021) need a corrected ",
                 "test (estimator = \"WLSMV\", lavaan's default for ",
                 "ordered data).")
        } else "and its standard measures are the ones to report."), arg = NULL)
    }
    if (version == "scaled" && !has_scaled) {
      .dmar_refuse_argument(c(
        "This fit has no scaled test statistic, so there are no scaled ",
        "measures; refit it with a robust estimator (for example estimator = ",
        "\"MLR\") or use version = \"standard\"."), arg = NULL)
    }
  }
  # The remedy the two warnings below name: the robust measures when the
  # fit has them, and otherwise a refit with a corrected test, since the
  # default "auto" can report only what lavaan computed. The standard
  # statistic of full WLS is chi square asymptotically, so its standard
  # measures, the only ones it has, need no warning. The standard
  # statistic of a ULS fit of continuous data has a warning of its own
  # below, since it is in the units of the squared covariances as well.
  if (version == "standard" && not_chisq && categorical) {
    .dmar_warn(c(
      "The standard statistic of a fit to ordered data is not chi square ",
      "distributed, so its p-value is NA and the RMSEA interval and the ",
      "close fit tests are not valid; ",
      if (has_robust) {
        "the default version = \"auto\" reports the corrected measures of Savalei (2021)."
      } else {
        paste0("refit with a corrected test (estimator = \"WLSMV\", lavaan's ",
               "default for ordered data) for the corrected measures of ",
               "Savalei (2021).")
      }), "result_qualified")
  } else if (version == "standard" && not_chisq && !uls_continuous) {
    .dmar_warn(c(
      "The standard statistic of a ", options$estimator, " fit is not ",
      "chi square distributed, so its p-value is NA and the RMSEA ",
      "interval and the close fit tests are not valid; ",
      if (has_robust) {
        "the default version = \"auto\" reports the robust measures."
      } else {
        paste0("leave out test = \"standard\", so that lavaan tests the ",
               "fit by Browne's residual based statistic, or refit with a ",
               "corrected test (estimator = \"WLSMV\" with ordered = ",
               "FALSE).")
      }), "result_qualified")
  }
  if (version == "scaled" && categorical) {
    .dmar_warn(c(
      "With ordered data the scaled measures do not estimate the ",
      "population RMSEA and CFI (Savalei, 2021); the default version = ",
      "\"auto\" reports the corrected ones."), "result_qualified")
  }
  if (version == "standard" && weighted) {
    .dmar_warn(c(
      "The standard statistic of a sampling weighted fit is a weighted ",
      "pseudo likelihood ratio statistic and is not chi square ",
      "distributed, so its p-value, RMSEA interval, and close fit tests ",
      "are not valid; ",
      if (has_robust) {
        "the default version = \"auto\" reports the robust measures."
      } else {
        paste0("refit with a corrected test (estimator = \"MLR\", or leave ",
               "out test = \"standard\" so that lavaan's default for ",
               "sampling weighted data applies) for the robust measures.")
      }), "result_qualified")
  }

  # The lavaan name behind each row under that version.
  sfx <- if (version == "standard") "" else paste0(".", version)
  test_sfx <- if (version == "standard" || !has_scaled) "" else ".scaled"
  key <- c(chi_square = paste0("chisq", test_sfx), df = paste0("df", test_sfx),
           p_chi_square = paste0("pvalue", test_sfx),
           rmsea = paste0("rmsea", sfx),
           rmsea_ci_lower = paste0("rmsea.ci.lower", sfx),
           rmsea_ci_upper = paste0("rmsea.ci.upper", sfx),
           p_close_fit = paste0("rmsea.pvalue", sfx),
           p_not_close_fit = paste0("rmsea.notclose.pvalue", sfx),
           cfi = paste0("cfi", sfx), tli = paste0("tli", sfx),
           nfi = paste0("nfi", sfx), ifi = paste0("ifi", sfx),
           rni = paste0("rni", sfx))
  if (version == "robust") key[c("nfi", "ifi")] <- NA_character_
  pick <- function(k) {
    if (!is.na(k) && k %in% names(fm)) .sem_number(fm[[k]]) else NA_real_
  }
  vals <- stats::setNames(rep(NA_real_, length(.SEM_FIT_TERMS)), .SEM_FIT_TERMS)
  for (term in names(key)) vals[[term]] <- pick(key[[term]])
  # After normal theory two-stage estimation the test of exact fit is the
  # residual based statistic of Savalei and Bentler (2009), asymptotically
  # chi square, as in fa_confirmatory(); lavaan's scaled statistic is an
  # approximation. After robust two-stage estimation it is lavaan's scaled
  # statistic (Savalei & Falk, 2014). The robust measures are lavaan's
  # either way, computed from the naive statistic and its expectation.
  if (identical(options$missing, "two.stage")) {
    rt <- .ts_residual_test(fit)
    vals[c("chi_square", "df", "p_chi_square")] <-
      c(rt$statistic, rt$df, rt$p_value)
  }
  # With auxiliary variables the SRMR is over the model's own variables
  # (.sem_srmr_own()), as the admissibility check sets them aside.
  srmr_aux <- .sem_aux_names(fit)
  vals[["srmr"]] <- if (length(srmr_aux)) {
    .sem_srmr_own(fit, srmr_aux)
  } else if ("srmr_bentler_nomean" %in% names(fm)) {
    pick("srmr_bentler_nomean")
  } else pick("srmr")
  vals[c("close_fit_bound", "not_close_fit_bound")] <- c(close_fit, not_close_fit)
  vals[c("N", "n_groups")] <- c(conv$N, conv$groups)
  df_used <- vals[["df"]]

  # The rows DMAR computes rather than reads. The scaled version is the
  # scaled statistic put in the standard formulas, its RMSEA, interval, and
  # close fit tests all from the chi_square and df rows (Brosseau-Liard et
  # al., 2012, Eqs. 11 and 13), as EQS and Mplus report them. lavaan 0.7-3
  # computes those only for a scaled and shifted statistic or under
  # rmsea.scaled.mplus = TRUE; otherwise it refers the unscaled statistic
  # to the noncentral chi square on tr(UGamma) = c df degrees of freedom,
  # which gives the same point estimate but another interval and other
  # p-values than the rows beside them imply.
  a <- conv$multiplier
  G <- conv$groups
  scale_c <- pick("chisq.scaling.factor")
  computable <- is.finite(df_used) && df_used > 0 && is.finite(vals[["chi_square"]])
  big <- computable &&
    max(vals[["chi_square"]],
        a * df_used * max(close_fit, not_close_fit)^2 / G) > .SEM_NCP_PCHISQ_LIMIT
  # The robust measures of the Brosseau-Liard family: lavaan's for a
  # Satorra-Bentler, Yuan-Bentler, or scaled and shifted statistic, and
  # DMAR's for the mean and variance adjusted one (`mva`).
  bl_family <- has_scaled && !fiml && !two_stage && !categorical &&
    (isTRUE(scaled_test %in% c("satorra.bentler", "yuan.bentler",
                               "yuan.bentler.mplus", "scaled.shifted",
                               "scaled.shifted.corrected")) ||
       length(mva) > 0L) &&
    is.finite(scale_c)
  rmsea_terms <- c("rmsea", "rmsea_ci_lower", "rmsea_ci_upper", "p_close_fit",
                   "p_not_close_fit")
  if (computable &&
      (version == "scaled" || (version == "standard" && big && !uls_units))) {
    vals[rmsea_terms] <- .sem_rmsea_rows(vals[["chi_square"]], df_used, a, G,
                                         conf_level, close_fit, not_close_fit)
  }
  # A statistic in the millions: lavaan's limits and close fit tests come
  # from R's pchisq(), which fails there (R/ci_nc_chisq_internals.R). The
  # robust ones of a complete data robust ML fit are recomputed from the
  # scaled statistic and its scaling constant (Brosseau-Liard et al., 2012,
  # Eq. 10, the interval lavaan forms); those of a mean and variance
  # adjusted statistic were computed so by .sem_robust_mva(); the others
  # cannot be, and are NA with a warning when the statistic is far enough
  # out to matter.
  if (version == "robust" && big && !length(mva)) {
    if (bl_family) {
      rows <- .sem_rmsea_rows(vals[["chi_square"]], df_used, a, G, conf_level,
                              close_fit, not_close_fit, c_hat = scale_c)
      vals[setdiff(rmsea_terms, "rmsea")] <- rows[setdiff(rmsea_terms, "rmsea")]
    } else {
      # The interval inverts a statistic of the size of the chi_square row,
      # and each test refers it to the noncentrality of its bound.
      ncp_bound <- a * df_used * c(close_fit, not_close_fit)^2 / G
      lost <- c(rmsea_ci_lower = vals[["chi_square"]] > 1e6,
                rmsea_ci_upper = vals[["chi_square"]] > 1e6,
                p_close_fit = ncp_bound[1L] > 1e6,
                p_not_close_fit = ncp_bound[2L] > 1e6)
      if (any(lost)) {
        vals[names(lost)[lost]] <- NA_real_
        what <- c(if (lost[["rmsea_ci_lower"]]) "RMSEA interval",
                  if (lost[["p_close_fit"]]) "test of close fit",
                  if (lost[["p_not_close_fit"]]) "test of not close fit")
        several <- length(what) > 1L
        .dmar_warn(c(
          "lavaan computes the robust ", paste(what, collapse = " and "),
          " of this fit from a noncentral chi square whose statistic or ",
          "noncentrality exceeds 10^6, where R's pchisq() loses its ",
          "accuracy, and DMAR cannot rebuild the corrected statistic of a ",
          "FIML, ordered data, or two-stage fit to recompute ",
          if (several) "them" else "it", ", so ",
          if (several) "they are" else "it is", " NA."), "result_qualified")
      }
    }
  }

  # A ULS fit of continuous data whose measures are formed from the ULS
  # statistic, (N - G) times the ULS fit function, in the units of the squared
  # covariances: every version but the standard one of a fit tested by
  # Browne's residual based statistic (estimator = "ULS" alone). A measure
  # is reported only where it does not change when every variable is
  # multiplied by the same number (.sem_uls_continuous() gives the reasons
  # and the sizes): the RMSEA rows of the robust and standard versions, and
  # the standard CFI, TLI, RNI, and IFI, which set the degrees of freedom
  # against the statistic, are NA; and when the model restricts the means,
  # which makes the estimates themselves change with the units, every RMSEA
  # row and incremental index of every version is NA. The robust CFI, TLI,
  # and RNI, the scaled measures, and the NFI of a fit that leaves the means
  # free are kept. The standard version of a fit tested by Browne's
  # statistic is reported where that statistic has no units: when the
  # estimates do not change with the units either, because the model
  # leaves the means free or restricts them through parameters the
  # covariances do not share (.sem_uls_estimates_units()); otherwise its
  # RMSEA rows and incremental indices are NA too (the scalar rung of the
  # first example of measurement_invariance() under estimator = "ULS" has a
  # CFI of .8874 for the raw scores and .8861 with every score divided by
  # 10). A saturated fit is left to the conventions
  # below, under which each index takes its value at a statistic of zero.
  uls_stat <- uls_continuous && (version != "standard" || uls_units) &&
    is.finite(df_used) && df_used > 0
  uls_browne <- uls_continuous && !uls_stat && version == "standard" &&
    is.finite(df_used) && df_used > 0
  uls_restricted <- if (uls_stat) {
    .sem_means_restricted(fit)
  } else if (uls_browne) {
    .sem_uls_estimates_units(fit)
  } else FALSE
  uls_case <- if (uls_browne) {
    if (isFALSE(uls_restricted)) "none" else "browne"
  } else if (!uls_stat) {
    "none"
  } else if (!isFALSE(uls_restricted)) {
    "means"
  } else if (version == "robust") {
    "robust"
  } else if (version == "standard") {
    "standard"
  } else "none"
  incremental <- c("cfi", "tli", "rni", "nfi", "ifi")
  uls_withheld <- switch(uls_case,
                         means = c(rmsea_terms, incremental),
                         browne = c(rmsea_terms, incremental),
                         robust = rmsea_terms,
                         standard = c(rmsea_terms, "cfi", "tli", "rni", "ifi"),
                         character(0))
  vals[uls_withheld] <- NA_real_
  uls_note <- switch(uls_case,
    means = paste0(
      if (is.na(uls_restricted)) {
        paste0("whether the model restricts the means could not be read, ",
               "and if it does, ")
      }, .sem_uls_means_text(), "; the RMSEA rows and the CFI, TLI, RNI, ",
      "NFI, and IFI are NA"),
    browne = paste0(
      .sem_uls_browne_text(uls_restricted), "; the RMSEA rows and the CFI, ",
      "TLI, RNI, NFI, and IFI are NA, and the chi_square row and the SRMR ",
      "describe this fit, whose estimates depend on the units"),
    robust = paste0(
      "lavaan's robust RMSEA of a ULS fit of continuous data changes with ",
      "the units of the variables, and the RMSEA rows are NA; the robust ",
      "CFI, TLI, and RNI do not change when every variable is multiplied ",
      "by the same number"),
    standard = if (isTRUE(options$correlation)) {
      paste0(
        "the standard statistic is (N - G) times the ULS discrepancy of the ",
        "correlation structure, which has no chi square reference ",
        "distribution, so the RMSEA rows and the CFI, TLI, RNI, and IFI, ",
        "which set the degrees of freedom against it and rest on the ML ",
        "discrepancy for their population meaning, are NA; the NFI, a ratio ",
        "of two such statistics, is reported")
    } else {
      paste0(
        "the standard statistic is (N - G) times the ULS fit function, in the ",
        "units of the squared covariances, so the RMSEA rows and the CFI, TLI, ",
        "RNI, and IFI, which set the degrees of freedom against it, change ",
        "with the units of the variables and are NA; the NFI, a ratio of two ",
        "such statistics, does not")
    },
    NA_character_)
  if (uls_continuous && version == "scaled" && uls_case == "none") {
    uls_note <- paste0("Scaled descriptive indices put the corrected ULS ",
                      "statistic into the index formulas. They do not estimate ",
                      "ML-population fit, and conventional ML cutoffs are not ",
                      "established for this use (Savalei, 2018).")
  }
  if (uls_case == "means") {
    .dmar_warn(c(
      "This ULS fit of continuous data ",
      if (isTRUE(uls_restricted)) "restricts the means, and " else {
        paste0("may restrict the means (DMAR could not read whether it ",
               "does), and if it does, ")
      },
      .sem_uls_means_text(), ". The RMSEA, its interval, the tests of close ",
      "and not close fit, and the CFI, TLI, RNI, NFI, and IFI are therefore ",
      "NA; the chi_square row and the SRMR, which describe this fit, are ",
      "reported. Maximum likelihood does not depend on the units (estimator ",
      "= \"MLR\" for nonnormal data)."), "result_qualified")
  } else if (uls_case == "browne") {
    .dmar_warn(c(
      .sem_uls_browne_text(uls_restricted), ". The RMSEA, its interval, the ",
      "tests of close and not close fit, and the CFI, TLI, RNI, NFI, and IFI ",
      "are therefore NA. The chi_square row, whose chi square reference ",
      "distribution holds at any consistent estimates (Browne, 1984, ",
      "Proposition 4), and the SRMR are reported; they describe this fit, whose ",
      "estimates depend on the units. Maximum likelihood does not depend on ",
      "the units (estimator = \"MLR\" for nonnormal data)."),
      "result_qualified")
  } else if (uls_case == "robust") {
    .dmar_refuse_argument(c(
      "ML-derived robust indices are unavailable for continuous ULS. ",
      "Use eligible scaled indices with their documented qualification."),
      arg = NULL)
  } else if (uls_case == "standard") {
    .dmar_warn(c(if (isTRUE(options$correlation)) {
              paste0("The standard statistic of a ULS fit of a correlation ",
                     "structure is (N - G) times the ULS discrepancy of the ",
                     "correlations: it is not chi square distributed (lavaan ",
                     "gives it no p-value), and the RMSEA, CFI, TLI, RNI, and ",
                     "IFI formed from it rest on the ML discrepancy for their ",
                     "population meaning, so they are NA, with the RMSEA ",
                     "interval and the tests of close and not close fit; the ",
                     "NFI, a ratio of two such statistics, is reported. ")
            } else {
              paste0("The standard statistic of a ULS fit of continuous data is (N - G) ",
                     "times the ULS fit function, which is in the units of the squared ",
                     "covariances: it is not chi square distributed (lavaan gives it no ",
                     "p-value), and the RMSEA, CFI, TLI, RNI, and IFI formed from it ",
                     "set its degrees of freedom, a number, against it, so they change ",
                     "with the units of the variables and are NA, with the RMSEA ",
                     "interval and the tests of close and not close fit; the NFI, a ",
                     "ratio of two such statistics, does not and is reported. ")
            },
            if (has_scaled) {
              paste0("The default version = \"auto\" reports eligible scaled ",
                     "indices from the corrected ULS test.")
            } else {
              paste0("Leave out test = \"standard\", so that lavaan tests the ",
                     "ULS fit by Browne's residual based statistic, or fit with ",
                     "estimator = \"ULSM\" or \"ULSMV\" for a corrected ",
                     "statistic.")
            }), "uls_standard_statistic")
  }

  # The boundary conventions: zero degrees of freedom, and no upper root.
  boundary <- character(0)
  exact <- TRUE
  if (is.finite(df_used) && df_used == 0) {
    vals[c("p_chi_square", "rmsea", "rmsea_ci_lower", "rmsea_ci_upper",
           "p_close_fit", "p_not_close_fit", "tli")] <- NA_real_
    # The statistic of a saturated model is zero in every version, and
    # with it the model term of the CFI, RNI, NFI, and IFI (T - c df, or
    # T), whatever the scaling constant c, which is 0/0 there: lavaan's
    # robust and scaled indices come out NA or not by rounding. Each takes
    # its value at T = 0, which lavaan's standard version computes (CFI 1,
    # by Bentler's, 1990, convention also when the baseline term is zero);
    # NFI and IFI stay NA in the robust version, which has none. A model
    # with no degrees of freedom whose statistic is not zero, as when an
    # inequality bound binds (.sem_exact_fit()), does not fit exactly: its
    # CFI, RNI, and IFI are those of the version, formed from its
    # statistic, and its standard NFI is 1 - T/T_b, which lavaan sets to 1
    # whenever df = 0 (.925, .925, .926, and .926 on one factor of three
    # variables correlated .8, .8, and .5 with bounds = "pos.var").
    exact <- .sem_exact_fit(fm)
    if (exact) {
      vals[c("rmsea", "rmsea_ci_lower", "rmsea_ci_upper")] <- 0
      vals[["tli"]] <- 1
      vals[["cfi"]] <- 1
      vals[["rni"]] <- pick("rni")
      if (version != "robust") {
        vals[c("nfi", "ifi")] <- c(pick("nfi"), pick("ifi"))
      }
    } else if (version == "standard") {
      vals[["nfi"]] <- .sem_ratio(pick("baseline.chisq") - pick("chisq"),
                                  pick("baseline.chisq"))
    }
    boundary <- .sem_boundary_text(if (exact) "zero_df" else "zero_df_misfit")
    .sem_warn_zero_df(exact)
  } else if (isTRUE(vals[["rmsea_ci_upper"]] == 0)) {
    boundary <- .sem_boundary_text("no_upper_root")
  }
  if (version == "robust" && is.finite(df_used) && !uls_continuous &&
      ((df_used > 0 &&
          (!is.finite(vals[["rmsea"]]) || !is.finite(vals[["cfi"]]))) ||
       (df_used == 0 && !exact && !is.finite(vals[["cfi"]])))) {
    if (df_used == 0) {
      .dmar_warn(c(
        "lavaan did not compute the robust CFI for this fit, which has ",
        "zero degrees of freedom and does not fit exactly (its scaling ",
        "constant is 0/0), so it is NA; version = \"standard\" gives ",
        "the standard CFI, formed from the statistic."), "result_qualified")
    } else {
      .dmar_warn(c(
        "lavaan did not compute the robust RMSEA or CFI for this fit ",
        "(with ordered data this happens when the polychoric ",
        "correlation matrix is not positive definite), so they are NA; ",
        "version = \"scaled\" gives the scaled versions, which do not ",
        "estimate the population values consistently."), "result_qualified")
    }
  }

  # Gamma hat: the standard version of a single group continuous ML fit of
  # complete, unweighted data whose means are not restricted, and the robust
  # version of such a fit with a Satorra-Bentler, Yuan-Bentler, scaled and
  # shifted, or mean and variance adjusted statistic, for which the robust
  # noncentrality is s (T_s - df_s) = T - c df for the statistic's scaling
  # constant s, its degrees of freedom df_s, and c = tr(U Gamma)/df
  # (lavaan's gfi_lrt.robust for the first three). The adjusted index and
  # the limits use the model's degrees of freedom, which differ from df_s
  # for the mean and variance adjusted statistic.
  gamma_ok <- version %in% c("standard", "robust") &&
    (version == "standard" || bl_family) &&
    identical(options$estimator, "ML") && G == 1L && !categorical && !fiml &&
    !two_stage && !weighted && !isTRUE(options$conditional.x)
  means_restricted <- if (gamma_ok) .sem_means_restricted(fit) else FALSE
  if (gamma_ok && isFALSE(means_restricted)) {
    excess <- if (version == "robust") {
      scale_c * (vals[["chi_square"]] - df_used) / a
    } else (vals[["chi_square"]] - df_used) / a
    gamma <- .sem_gamma_hat(excess, pick("df"), length(lavaan::lavNames(fit, "ov")),
                            vals[["rmsea_ci_lower"]], vals[["rmsea_ci_upper"]])
    vals[names(gamma)] <- gamma
  }

  # Cases removed by listwise deletion. (Under sampling weights lavaan's N
  # is the sum of the rescaled weights, equal to the number of cases up to
  # rounding error, hence the comparison at one half.)
  # (Under FIML lavaan drops only the cases that observe no variable, which
  # carry no information, so no warning. Under two-stage estimation it
  # does the same, except that a case observed on an auxiliary variable
  # alone is read by its first stage and still left out of its N, which
  # the second warning counts.)
  n_orig <- sum(lavaan::lavInspect(fit, "norig"))
  if (!fiml && !two_stage && n_orig - conv$N > 0.5) {
    .dmar_warn(
      sprintf(paste0("lavaan removed %d of %d cases by listwise deletion; ",
                     "the measures use the %d complete cases, which assumes ",
                     "the data are missing completely at random. Refit with ",
                     "missing = \"ml\" to use every case."),
              round(n_orig - conv$N), round(n_orig), round(conv$N)),
      "result_qualified")
  }
  # The second warning (.ts_aux_only_warning() in R/two_stage_internals.R,
  # which sem_information_criteria() and compare_sem() give as well).
  .ts_aux_only_warning(
    fit,
    rest = paste0("the naive statistic, and with it lavaan's scaled ",
                  "statistic, the corrected indices, the RMSEA multiplier, ",
                  "and the N row, rest"),
    consequence = "the rows do not describe one set of cases",
    aside = if (identical(options$missing, "two.stage")) {
      "; the chi_square row, the residual based statistic, uses every case"
    } else "")

  test <- if (two_stage) {
    .ts_test_label(options$missing)
  } else if (version != "standard" && has_scaled) {
    paste0(scaled_test, " (scaled)")
  } else standard_test
  correction <- if (version == "robust") {
    if (two_stage) {
      paste0("two-stage corrected indices of Zhang and Savalei (2023), ",
             if (identical(options$missing, "robust.two.stage")) "TS V2"
             else "the normal theory analog of TS V2")
    } else if (fiml) {
      paste0("FIML corrected indices of Zhang and Savalei (2023), version ",
             if (has_scaled) "3 (no normality assumption)" else "6 (normal theory)")
    } else if (categorical) {
      "categorical correction of Savalei (2021)"
    } else if (isTRUE(scaled_test %in% c("scaled.shifted", "scaled.shifted.corrected"))) {
      "mean and variance adjusted correction of Savalei (2018)"
    } else if (length(mva) && !is.null(attr(mva, "scaling"))) {
      # Only where DMAR computed them from the constants: a saturated fit
      # gets the values at a statistic of zero (.sem_robust_mva()), with
      # nothing computed from c, which is 0/0 there, and the plain wording
      # below, which a saturated MLM fit gets too.
      paste0("nonnormality corrections of Brosseau-Liard, Savalei, and Li (2012) ",
             "and Brosseau-Liard and Savalei (2014), computed by DMAR from the ",
             "scaling constant and degrees of freedom of the mean and variance ",
             "adjusted statistic, for which lavaan computes none (Savalei, 2018)")
    } else if (uls_continuous) {
      paste0("lavaan's nonnormality corrections of Brosseau-Liard, Savalei, and ",
             "Li (2012) and Brosseau-Liard and Savalei (2014), derived for the ML ",
             "fit function, applied to the ULS statistic (Savalei, 2018)")
    } else {
      paste0("nonnormality corrections of Brosseau-Liard, Savalei, and Li (2012) ",
             "and Brosseau-Liard and Savalei (2014)")
    }
  } else NA_character_
  details <- list(
    source = "fitted lavaan model",
    version_requested = version_requested, version = version,
    version_reason = reason, test = test, correction = correction,
    estimator = options$estimator, likelihood = options$likelihood,
    missing = options$missing,
    two_stage_rows = ts_rows[c("read", "nobs", "extra", "corrected", "reestimated")],
    sample_size_convention = conv$convention, multiplier = conv$multiplier,
    N = conv$N, groups = conv$groups,
    srmr = if (length(srmr_aux)) {
      paste0("covariance residuals, Bentler standardization, over the model's ",
             "own variables (auxiliary variables set aside: ",
             paste(srmr_aux, collapse = ", "), ")")
    } else "covariance residuals, Bentler standardization (lavaan's srmr_bentler_nomean)",
    gamma = if (isTRUE(means_restricted)) {
      paste0("not computed: the model restricts the means, and gamma hat is ",
             "defined for covariance structures")
    } else if (is.na(means_restricted)) {
      paste0("not computed: lavaan's derivatives of the moments could not be ",
             "read to tell whether the model restricts the means")
    } else if (gamma_ok) "computed" else {
      paste0("not computed: defined here for the standard version, and the robust ",
             "version of a robust ML fit of complete data, of a single group ",
             "continuous ML fit without FIML or sampling weights")
    },
    uls = uls_note,
    boundary = boundary)
  details$footer <- c(
    .ts_rows_line(ts_rows, kind = "cfa", lavaan_missing = options$missing),
    paste0("Version: ", version, " (",
           if (version_requested == "auto") paste0("auto: ", reason) else reason,
           ")"),
    paste0("Chi square: ", test,
           if (!is.na(correction)) paste0("; robust measures: ", correction)),
    if (!is.na(uls_note)) paste0("ULS of continuous data: ", uls_note),
    sprintf("Sample size: N = %s in %s; multiplier a = %s (%s)",
            format(conv$N),
            if (conv$groups == 1L) "1 group" else paste(conv$groups, "groups"),
            format(conv$multiplier),
            if (conv$convention == "N") "lavaan's normal likelihood forms T with N"
            else "lavaan forms T with N less the number of groups"))
  list(values = vals, details = details)
}

.sem_boundary_text <- function(which) {
  switch(which,
    zero_df = .sem_saturated_text(),
    zero_df_misfit = paste0(
      "zero degrees of freedom and a statistic that is not zero (an ",
      "inequality constraint that binds costs no degree of freedom in ",
      "lavaan's count), so the model does not fit exactly: the exact fit ",
      "p-value, the RMSEA and its limits, the close fit tests, and TLI are ",
      "undefined and NA, and the CFI, RNI, NFI, and IFI are formed from ",
      "the statistic in the versions that have them"),
    no_upper_root = paste0("the statistic lies at or below the lower ",
                           "(1 - conf_level)/2 quantile of the central chi ",
                           "square, so no noncentrality is consistent with it ",
                           "and the upper RMSEA limit is 0"))
}

.sem_warn_zero_df <- function(exact = TRUE) {
  if (exact) {
    .dmar_warn(.sem_saturated_text(), "result_qualified")
  } else {
    .dmar_warn(c(
      "The model has zero degrees of freedom but its statistic is not ",
      "zero, so it does not fit exactly (an inequality constraint that ",
      "binds, such as a variance held at a bound, costs no degree of ",
      "freedom in lavaan's count, and a model that is not identified ",
      "at its estimates can misfit too). The exact fit p-value, the ",
      "RMSEA and its limits, the tests of close and not close fit, and ",
      "TLI divide by the degrees of freedom and are undefined and ",
      "returned as NA; the CFI, RNI, NFI, and IFI are formed from the ",
      "statistic in the versions that have them."), "result_qualified")
  }
}
