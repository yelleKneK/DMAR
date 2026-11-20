#' Measurement Invariance Across Groups
#'
#' Fits the standard ladder of multiple group invariance models for a
#' measurement model of any number of factors and reports the comparison
#' of its rungs in one table: configural invariance (same pattern, all
#' parameters free by group), metric (equal loadings; required before
#' comparing relations involving the factor), scalar (equal loadings and
#' intercepts; required before comparing factor means), and strict (equal
#' residual variances as well; required before comparing observed-score
#' variances). With ordered indicators a thresholds rung comes first,
#' because thresholds rather than intercepts carry the location information
#' (Wu & Estabrook, 2016); with dichotomous indicators the two are not
#' separately identified, so that rung is folded into metric. Each rung is
#' tested against the previous with a likelihood ratio test, scaled when the
#' estimator in force carries a robust test, which is what declaring ordered
#' items arranges, and the practical-fit changes (delta CFI, delta RMSEA)
#' are reported alongside, since with large samples the chi square will flag
#' trivial differences (Cheung & Rensvold, 2002, suggest delta CFI of about
#' -.01 as a red flag). The fitted \pkg{lavaan} objects come back as an
#' attribute, so score tests and partial invariance refits do not require
#' refitting the ladder. Requires \pkg{lavaan}.
#'
#' @param data A \code{data.frame} with the items and the grouping
#'   variable.
#' @param model The measurement model, given either as \pkg{lavaan} model
#'   syntax (a single string, or a character vector of lines that is
#'   collapsed with newlines, such as
#'   \code{c("visual =~ x1 + x2 + x3", "verbal =~ x4 + x5 + x6")}) or as a
#'   named list mapping each factor name to a character vector of its item
#'   names, from which the syntax is built. Exactly one of \code{model} and
#'   \code{items} must be supplied. An indicator of either form that
#'   \code{data} lacks is refused with an error of class
#'   \code{dmar_items_invalid} (see \code{\link{dmar_metadata}}).
#' @param group Single character string naming the grouping column (two or
#'   more groups).
#' @param items Character vector (three or more) naming the indicator
#'   columns of a single factor. A convenience for the one-factor
#'   (congeneric) case of \code{\link{cfa_1}}, equivalent to
#'   \code{model = "f =~ item_1 + item_2 + ..."}. Exactly one of
#'   \code{model} and \code{items} must be supplied. Items that cannot be
#'   read from \code{data} (a name that is not a column, a name given
#'   twice, fewer than three) are refused with an error of class
#'   \code{dmar_items_invalid}.
#' @param levels Which rungs of the ladder to fit, in order: any leading
#'   subset of the ladder in force. With continuous indicators the ladder is
#'   \code{c("configural", "metric", "scalar", "strict")}; with ordered
#'   indicators it is \code{c("configural", "thresholds", "metric",
#'   "scalar", "strict")}; when every indicator is ordered and dichotomous it
#'   is \code{c("configural", "metric", "scalar")}, because neither the
#'   thresholds nor the strict rung is identified for a two-category item
#'   (see \code{ordered} below). \code{NULL} (the default) fits the whole
#'   ladder in force.
#' @param ordered Ordered categorical (including binary) items: \code{TRUE}
#'   for every indicator in the model, or a character vector naming the
#'   ordered ones. \code{NULL} (default) treats all indicators as
#'   continuous. Declaring ordered items switches the ladder and, unless the
#'   estimator already carries the mean and variance adjusted test, switches
#'   the estimator (see \code{estimator} below). The number of observed
#'   categories per declared item is counted before the ladder is built,
#'   because a dichotomous item has a single threshold that is not
#'   separately identified from the intercept of its underlying response
#'   (Millsap & Yun-Tein, 2004; Wu & Estabrook, 2016), and its residual
#'   variance is not a free parameter under either parameterization. When
#'   every declared item is dichotomous the thresholds rung is folded into
#'   metric, and the strict rung is dropped as well unless a continuous
#'   indicator is left to carry a residual variance; when only some declared
#'   items are dichotomous, both rungs are kept and a message names those
#'   items, whose thresholds and residual variances the rungs leave untested.
#' @param estimator Estimator passed to \pkg{lavaan}. Defaults to
#'   \code{"ML"}; use \code{"MLR"} for robust corrections, in which case the
#'   likelihood ratio tests use the scaled difference. When \code{ordered}
#'   is supplied, a maximum likelihood or generalized least squares (the \code{"GLS"} label; Browne, 1974)
#'   estimator is replaced by \code{"WLSMV"}, and \code{"DWLS"} and
#'   \code{"ULS"} are replaced by \code{"WLSMV"} and \code{"ULSMV"}. That
#'   second substitution changes only the test: in \pkg{lavaan}
#'   \code{"WLSMV"} is \code{"DWLS"} and \code{"ULSMV"} is \code{"ULS"} with
#'   the mean and variance adjusted statistic, which is what makes the
#'   rung-to-rung difference test asymptotically valid for ordered data. The
#'   one categorical estimator left alone is \code{"WLS"}, the full weight
#'   matrix estimator of Browne (1984) and Muthén (1984), whose statistic is
#'   asymptotically chi square already and whose weight matrix a
#'   substitution would silently change.
#' @param missing Missing data handling passed to \pkg{lavaan}: one of
#'   \code{"fiml"} (the default) or \code{"ml"} (full information maximum
#'   likelihood, which fits every case with an observed item),
#'   \code{"listwise"} (listwise deletion), \code{"pairwise"},
#'   \code{"two_stage"}, or \code{"robust_two_stage"} (two-stage maximum
#'   likelihood, with normal theory or nonnormality robust standard errors
#'   and tests; lavaan's \code{"two.stage"} and \code{"robust.two.stage"}
#'   are accepted for these two). Full information maximum likelihood and
#'   two-stage estimation are not available with the categorical
#'   estimator, and two-stage estimation needs \code{estimator = "ML"};
#'   asking for another combination is an error. The default analysis is
#'   full information maximum likelihood whenever a value of a model
#'   variable is missing and the complete data analysis otherwise, unlike
#'   \code{\link{fa_confirmatory}}, which fits complete data by FIML too:
#'   data in which no such value is missing, or none but in cases with
#'   every model variable missing (which both treatments set aside), are
#'   fitted as they stand (under FIML the robust indices of an MLR ladder
#'   would be the corrected ones of Zhang and Savalei, 2023, not the
#'   complete data ones), and incomplete data by FIML. FIML needs
#'   continuous items,
#'   \code{estimator = "ML"} or \code{"MLR"}, standard errors other than
#'   \code{se = "robust.sem"}, and a test other than the Satorra-Bentler
#'   kind passed through \code{...}; when the default meets incomplete
#'   items and one of these, the call stops and names the option to choose
#'   (\code{missing = "listwise"}, \code{"pairwise"} for ordered items,
#'   \code{estimator = "MLR"}, or, with \code{estimator = "MLR"}, leaving
#'   out \code{test}) rather than deleting cases unasked, and an explicit
#'   \code{missing = "fiml"} with \code{se = "robust.sem"} or such a test
#'   stops too, whether or not a value is missing, since \pkg{lavaan}
#'   would fit the complete cases instead. Under FIML a group that never
#'   observes a variable of the model, or never observes two of them
#'   together, is refused before any fit (an error of class
#'   \code{dmar_group_coverage}). See \emph{Details}.
#' @param group_partial Character vector of parameters to leave free across
#'   groups at every rung, passed to \pkg{lavaan}'s \code{group.partial}.
#'   Either user-supplied parameter labels or \pkg{lavaan} parameter
#'   specifications such as \code{"visual =~ x2"} or \code{"x2 ~1"}. This is
#'   the partial invariance case of Byrne, Shavelson, and Muthén (1989).
#'   \code{NULL} (default) constrains everything the rung calls for, as does
#'   an empty vector (\code{character(0)} or \pkg{lavaan}'s own default
#'   \code{""}), which frees nothing in \pkg{lavaan} either.
#' @param parameterization Identification of ordered indicators, passed to
#'   \pkg{lavaan}: \code{"delta"} (the default, lavaan's) or \code{"theta"}.
#'   Residual variances are free parameters only under \code{"theta"}, so a
#'   strict rung requested with ordered items switches to \code{"theta"}
#'   with a message. Ignored when no item is ordered.
#' @param ... Further arguments passed to \code{\link[lavaan]{cfa}} at every
#'   rung (for example \code{std.lv}, \code{orthogonal}, \code{cluster}).
#'   \pkg{lavaan}'s own \code{group.partial} is merged into
#'   \code{group_partial}, and an empty one (\code{NULL},
#'   \code{character(0)}, or \code{""}) frees nothing, as in \pkg{lavaan},
#'   so a call such as \code{group.partial = if (free) terms else NULL}
#'   runs either way; \code{group.equal}, which the ladder sets rung
#'   by rung, is an error here (choose the rungs with \code{levels}).
#'
#' @details
#' The models are nested by construction, each adding equality constraints
#' across groups to the previous. With continuous indicators the constraint
#' sets are none (beyond the configuration), then \code{group.equal =
#' "loadings"}, then \code{c("loadings", "intercepts")}, then
#' \code{c("loadings", "intercepts", "residuals")}.
#'
#' The ordered ladder differs, and the difference is substantive rather than
#' cosmetic. For an ordered indicator the observed response is a coarsening
#' of an underlying continuous response at a set of thresholds, and it is
#' the thresholds, not an intercept, that locate the item on the latent
#' scale. Constraining loadings while leaving thresholds free across groups
#' therefore does not deliver what metric invariance is supposed to deliver,
#' and the accepted sequence (Millsap & Yun-Tein, 2004; Wu & Estabrook,
#' 2016) constrains thresholds first: \code{"thresholds"}, then
#' \code{c("thresholds", "loadings")} for metric, then
#' \code{c("thresholds", "loadings", "intercepts")} for scalar, then adding
#' \code{"residuals"} for strict. The intercept rung is not vacuous even
#' when every indicator is ordered: once thresholds and loadings are
#' constrained, \pkg{lavaan} frees the underlying-response intercepts in the
#' non-reference groups, and the scalar rung is what returns them to zero
#' and lets the latent means be estimated instead. Residual variances,
#' however, are free parameters only under the theta parameterization, so
#' \code{parameterization} switches to \code{"theta"} when a strict rung is
#' requested with ordered items.
#'
#' Dichotomous items are the exception the ordered ladder has to make room
#' for. A two-category item contributes one threshold, and that threshold,
#' the intercept of the underlying response, and its residual variance are
#' not separately identified: the data give one proportion per group per
#' item, which pins down a single standardized location and nothing else
#' (Millsap & Yun-Tein, 2004; Wu & Estabrook, 2016). Constraining thresholds
#' alone across groups is then a reparameterization rather than a
#' restriction, since \pkg{lavaan} frees the underlying-response intercepts
#' by exactly as many parameters as the constraint removes, and the residual
#' variances stay fixed at one in every group whichever parameterization is
#' in force. So the function counts the observed categories of each declared
#' ordered indicator before building the ladder, and a message explains what
#' the count implies. When all of them are dichotomous the thresholds rung is
#' folded into metric, which constrains thresholds and loadings together; the
#' strict rung goes too when every indicator in the model was declared
#' ordered, leaving configural, metric, scalar. A continuous indicator
#' alongside dichotomous ones keeps the strict rung, since its residual
#' variance is a free parameter. When only some declared items are
#' dichotomous the full ordered ladder is kept, since the polytomous items
#' still carry testable thresholds and residual variances, and the message
#' names the dichotomous items so their contribution to those two rungs is
#' not overread. Every rung that is dropped is one that would have cost zero
#' degrees of freedom.
#'
#' When the estimator carries a robust test, the difference between two chi
#' square statistics is not itself chi square distributed.
#' \code{\link[lavaan]{lavTestLRT}} then returns the scaled difference test
#' (the Satorra-Bentler or Satorra correction, chosen by \pkg{lavaan} to
#' match the test in force), and that is what \code{delta_chi_square} and
#' \code{p_value} report; the \code{"test"} attribute names the test used. In
#' that case the model-level \code{chi_square}, \code{p_chi_square},
#' \code{cfi}, and \code{rmsea} columns are lavaan's scaled and robust
#' versions, so \code{delta_chi_square} will not equal the difference of
#' consecutive \code{chi_square} entries. That is a property of scaled
#' difference testing, not an inconsistency.
#'
#' One version of the test and one of the indices hold every rung, so that
#' \code{delta_cfi} and \code{delta_rmsea} compare like with like, and the
#' \code{"model_test"} and \code{"fit_indices"} attributes name the
#' versions the columns hold. Where \pkg{lavaan} does not compute the
#' robust indices at every rung, every rung falls back to one other
#' version, and a line under the table says which. For a maximum
#' likelihood fit of continuous items, whose robust indices estimate the
#' population values of the ML indices with their small sample bias
#' corrected for nonnormality, it is the standard indices, which estimate
#' the same population values with the bias corrected as if the data were
#' normal; the scaled indices, which put the scaled statistic into the
#' standard formulas, estimate population values that change with the
#' kurtosis of the data (Savalei, 2018). Beside the older mean and variance
#' adjusted statistic (\code{test = "mean.var.adjusted"}, passed through
#' \code{...}) \pkg{lavaan} computes no robust indices, and they are
#' computed here, as \code{\link{sem_fit_measures}} computes them: they do
#' not depend on which corrected statistic accompanies them (Savalei,
#' 2018), so they equal those of the Satorra-Bentler statistic of the same
#' fit, and a line under the table says that DMAR computed them. For
#' ordered items, whose standard statistic is not the basis of an index, it
#' is the scaled indices: \pkg{lavaan} leaves the robust indices \code{NA}
#' at the rungs that constrain the thresholds under the default delta
#' parameterization, as in the ordered example below. A ladder of
#' continuous items fit by ULS (\code{estimator = "ULS"}, \code{"ULSM"}, or
#' \code{"ULSMV"}) withholds the ML-derived robust RMSEA and incremental
#' indices (Savalei, 2018). Invariance to a common change of units alone
#' does not establish validity outside their derivation. Eligible scaled
#' indices are reported with a qualification: they do not estimate the
#' same population discrepancy as ML-based indices, and ML cutoffs do not
#' automatically apply. The standard ULS statistic is in the units of
#' squared covariances, so its CFI and RMSEA are withheld too. At a rung
#' that restricts the means, as the scalar and
#' strict rungs do through the intercepts, the ULS fit function adds the
#' squared mean residuals, in the squared units of the items, to the
#' squared covariance residuals, in the squared units of the covariances,
#' so the statistic, the loadings, and every index change with the units:
#' on the first example under \code{"ULSM"}, the scalar CFI is .9885 for the
#' raw scores and .9840 with every score divided by 10, and the scalar
#' delta CFI -.00005 and -.0046, a judgment against the cutoff of Cheung and
#' Rensvold (2002) that turns on the units. The \code{cfi} and \code{rmsea}
#' columns of such a rung are \code{NA}, and the \code{chi_square} column
#' and the difference test are kept, since their reference distribution
#' does not depend on the units, though their values do. A ladder fit by
#' \code{estimator = "ULS"} alone is tested by Browne's residual based
#' statistic, which has no units, and its \code{cfi} and \code{rmsea}
#' columns, formed from that statistic, are reported at the rungs whose ULS
#' estimates do not change with the units, those that leave the means free.
#' At a rung that restricts the means through parameters the covariances
#' share, as the scalar and strict rungs do through the loadings, which
#' enter the means through the latent means, the ULS estimates change with
#' the units, and so do Browne's statistic evaluated at them and the
#' indices formed from it: on the first example, the scalar CFI is .8874
#' for the raw scores and .8861 with every score divided by 10, and the
#' scalar delta CFI -.0212 and -.0225. Those columns are \code{NA} there
#' too, and the \code{chi_square} column is kept. \pkg{lavaan}'s
#' difference test of a ULS ladder without a corrected statistic refers the
#' difference of the ULS statistics, in the units of the squared
#' covariances and not chi square distributed, to the chi square
#' distribution, so \code{delta_chi_square} and \code{p_value} are
#' \code{NA} at every rung, and \code{"test"} is \code{NA}; \code{"ULSM"} and
#' \code{"ULSMV"} give a scaled difference test. \code{estimator = "MLR"},
#' which does not depend on the units, is the estimator for nonnormal
#' continuous items. When \pkg{lavaan}
#' reports no scaled statistic at every rung (beside a bootstrap test,
#' \code{test = "bollen.stine"}), the \code{chi_square} column and the
#' difference test are the standard ones, and the attributes say so. A
#' saturated rung (zero degrees of freedom and a statistic of zero) fits
#' exactly in every version
#' and does not decide, since \pkg{lavaan} can leave its robust CFI
#' \code{NA} there from a ratio of zeros, or not compute it. Its statistic
#' is zero in every version, and with it the model term of the CFI,
#' \eqn{T - c \cdot df},
#' whatever the scaling constant \eqn{c} (which is 0/0 there), so its
#' \code{cfi} is 1 in every version; its RMSEA, \eqn{\sqrt{(T - c \cdot
#' df)/(N \cdot df)}}, is 0/0, which no version defines, and it keeps
#' \pkg{lavaan}'s 0, the convention \pkg{lavaan} applies in every version
#' to a model that fits exactly, so that \code{delta_cfi} and
#' \code{delta_rmsea} at the next rung compare with exact fit. A rung with
#' zero degrees of freedom whose
#' statistic is not zero, as when an inequality bound binds (\pkg{lavaan}'s
#' count of the degrees of freedom leaves the bound out), does not fit
#' exactly: its \code{cfi} is formed from its statistic and decides the
#' version of the indices (a ladder whose robust CFI \pkg{lavaan} leaves
#' \code{NA} at such a rung falls back to another version at every rung,
#' the standard indices for a maximum likelihood ladder of continuous
#' items), while its \code{rmsea}, which divides by the degrees of freedom,
#' is undefined and \code{NA}, as is \code{delta_rmsea} at the next rung,
#' with a line under the table that says so; it does not decide the
#' version of the test, which no version defines on zero degrees of
#' freedom.
#'
#' Which estimators carry such a test is worth being precise about, because
#' a naive chi square difference on ordered data is not asymptotically valid.
#' \code{"MLM"}, \code{"MLMV"}, \code{"MLR"}, \code{"WLSM"}, \code{"WLSMV"},
#' \code{"ULSM"}, and \code{"ULSMV"} carry one. \code{"DWLS"} and
#' \code{"ULS"} do not, which is why declaring ordered items promotes them to
#' \code{"WLSMV"} and \code{"ULSMV"}: the discrepancy function and hence the
#' parameter estimates are untouched, and only the statistic changes. The
#' remaining case is \code{"WLS"}, the full weight matrix estimator, whose
#' statistic is asymptotically chi square under the theory of Browne (1984)
#' and Muthén (1984) and so needs no correction; the \code{"test"} attribute
#' reports the standard difference test there. The sample size that theory
#' asks for is large, which limits full weighted least squares to large
#' samples, but that is a separate matter from whether the
#' difference test is the right one.
#'
#' \strong{Full information maximum likelihood.} With \code{missing =
#' "fiml"} (the default when an item value is missing) or \code{"ml"}, so
#' that the ladder uses every case with an observed item, as
#' \code{\link{fa_confirmatory}} fits the same model by default, the
#' \code{chi_square} column is the FIML
#' likelihood ratio statistic (its scaled version with \code{estimator =
#' "MLR"}), and the \code{cfi} and \code{rmsea} columns are \pkg{lavaan}'s
#' robust versions, the corrected indices of Zhang and Savalei (2023) that
#' \code{\link{fa_confirmatory}} reports for the same model: the complete data fit
#' function evaluated at the FIML estimates, with the saturated FIML
#' estimates in place of the sample moments, so that they estimate the
#' values the indices would have had with complete data (\pkg{lavaan} uses
#' their version 6, which assumes normality, and their version 3, which
#' does not, with \code{estimator = "MLR"}). The usual indices computed
#' from the FIML statistic move with the amount and kind of missing data
#' (Zhang & Savalei, 2023), so their change from rung to rung is not the
#' change complete data would show, and a cutoff for a change in fit, such
#' as the delta CFI of about -.01 of Cheung and Rensvold (2002), is a
#' statement about the latter. \code{delta_cfi} and \code{delta_rmsea}
#' therefore compare the corrected indices, and one version holds every
#' rung: when \pkg{lavaan} cannot compute them at some rung with degrees
#' of freedom, every rung reports its scaled indices (for a fit with a
#' scaled statistic) or its standard ones, and the line under the table
#' says which. \pkg{lavaan} cannot combine FIML with
#' the Satorra-Bentler or the mean and variance adjusted statistic
#' (\code{test = "satorra.bentler"} or \code{"mean.var.adjusted"} passed
#' through \code{...}), nor with \code{se = "robust.sem"}: it fits the
#' complete cases instead, with a warning. Under the default of
#' \code{missing} such a ladder is fitted as it stands when no item value
#' is missing and stops otherwise; asked for by name, FIML with one of
#' them stops before any fit whether or not a value is missing, as in
#' \code{\link{fa_confirmatory}}, and so does a ladder whose fits
#' \pkg{lavaan} made by another treatment than the FIML asked for (it
#' fits \code{estimator = "GLS"}, \code{"WLS"}, and \code{"ULS"} in two
#' stages).
#'
#' FIML estimates the means, variances, and covariances of the model's
#' variables within each group, so every group must observe every
#' variable, and every pair of them together, at least once. A group
#' that does not is refused by name before any fit, in an error of class
#' \code{dmar_group_coverage} that carries the variables and pairs at
#' issue (see \code{\link{dmar_metadata}}); \pkg{lavaan} would stop with
#' an internal error. \code{holzinger_swineford} is such a data set, by
#' school: \code{t25_paper_form_board_r} and \code{t26_flags} were given
#' in Grant-White only. The message names what to do: leave the variable
#' (or one of the pair) out of the model, or an auxiliary variable out of
#' \code{aux}, for this comparison, so that the ladder compares the groups
#' on what every group observes, or, when at least two groups would
#' remain, leave the groups at issue out of \code{data}. Listwise
#' deletion is no remedy, since a group that never observes a variable,
#' or never observes two together, has no case observed on every
#' variable, and the listwise ladder would stop on that group instead.
#'
#' Auxiliary variables passed through \code{...} as \pkg{lavaan}'s
#' \code{aux} enter a FIML fit as saturated correlates (Graham, 2003),
#' and \pkg{lavaan} makes the mean and variance of each one ordinary
#' parameters of each group, which the intercepts and residuals rungs
#' would equate across the groups; the rung tests would then test whether
#' the auxiliary is distributed alike in the groups, which says nothing
#' about the measurement model. They are left free at every rung, beside
#' any \code{group_partial} terms, so the ladder has the degrees of freedom
#' of the ladder without the auxiliary. On six Holzinger and Swineford
#' tests by school with 60 sentence completion scores deleted at random
#' and general information (whose means are 44.8 and 36.6 in the two
#' schools) as the auxiliary, the scalar rung's difference test is 23.88
#' on 4 df and its delta CFI -.031 (24.03 and -.031 without the
#' auxiliary).
#'
#' \strong{Two-stage estimation.} With \code{missing = "two_stage"} every
#' rung is fit in two stages (Yuan & Bentler, 2000; Savalei & Bentler,
#' 2009): the means and covariances are first estimated in each group by
#' maximum likelihood under the saturated normal model (the EM estimates,
#' which use every partially observed case and any auxiliary variables
#' passed through \code{...} as lavaan's \code{aux}), and the rung is then
#' fit to them by the complete data maximum likelihood discrepancy
#' function. A case observed on an auxiliary variable but on none of the
#' model's variables is read by the first stage and kept, as in
#' \code{\link{fa_confirmatory}}, which says why: \pkg{lavaan}'s two-stage
#' estimator leaves such a case out of its sample size, so every rung is
#' computed again with its sample size in each group set to the cases the
#' first stage read there. The statistics, the indices, the difference
#' tests, and the fits of the \code{"fits"} attribute then count every
#' case. The second stage weights each group's discrepancy by its share of
#' the sample size, so a rung that constrains parameters across groups has
#' estimates that depend on the groups' counts, and when those grow in
#' different ratios the rung is estimated again; the configural rung, which
#' constrains nothing across groups, keeps \pkg{lavaan}'s estimates. On six
#' Holzinger and Swineford tests by school with the whole block missing
#' for 77 pupils at random given general information (the auxiliary), 130
#' and 94 cases with an observed test grow to 156 and 145 read; the common
#' intercept of lozenges at the scalar rung is 18.00 (SE .697) rather than
#' \pkg{lavaan}'s 18.21 (SE .782), and the scalar rung's difference test is
#' 14.97 on 4 df, against 10.86 from \pkg{lavaan}'s own ladder, 16.46
#' under FIML with the same auxiliary, and 17.44 after listwise deletion.
#' A message says how many such cases there are, and a line under the
#' table names the sample size. \pkg{lavaan}'s \code{summary()} of a fit
#' in the \code{"fits"} attribute counts them as used, and
#' \code{\link[lavaan]{lavTestLRT}} warns that a difference test between
#' such a fit and one \pkg{lavaan} made of the same data on its own sample
#' size is invalid. The standard errors are the two-stage sandwich, and the
#' \code{chi_square} column of each rung is
#' the residual based statistic of Savalei and Bentler (2009, Eq. 8),
#' whose rejection rate their simulations with incomplete normal data
#' found close to the nominal level, and the statistic
#' \code{\link{fa_confirmatory}} reports after two-stage estimation. With
#' \code{missing = "robust_two_stage"} the standard errors are robust to
#' nonnormality (Yuan & Bentler, 2000) and \code{chi_square} is lavaan's
#' scaled statistic (Savalei & Falk, 2014). Under both, \code{cfi} and
#' \code{rmsea} are lavaan's robust indices, the corrected indices of
#' Zhang and Savalei (2023): their TS V2 after robust two-stage estimation,
#' and its normal theory analog after normal theory two-stage estimation
#' (\code{\link{fa_confirmatory}} explains the difference), so \code{delta_cfi}
#' and \code{delta_rmsea} compare those, and the
#' rung-to-rung test is the scaled difference test of Satorra and Bentler
#' (2001), which \code{\link[lavaan]{lavTestLRT}} applies to two-stage fits.
#' The difference is formed from the two rungs' naive statistics and their
#' scaling, so \code{delta_chi_square} is not the difference of the
#' \code{chi_square} entries. Two-stage estimation needs continuous items
#' and \code{estimator = "ML"}; robustness to nonnormality comes from
#' \code{"robust_two_stage"}, not from \code{estimator = "MLR"}. The
#' residual based statistic treats every observed variable as random, so a
#' model given as syntax that regresses on an observed variable needs
#' \code{fixed.x = FALSE} passed through \code{...}: under
#' \code{"two_stage"} the statistic is refused for a fit that conditions
#' on such a variable, as it is for sampling weights, clustered cases,
#' inequality constraints, and a parameter on a bound, each with a message
#' that says what to do instead.
#'
#' A rung whose constraints cost no degrees of freedom tests nothing, so
#' \code{delta_chi_square} and \code{p_value} are \code{NA} when
#' \code{delta_df} is zero rather than reporting a difference in chi square
#' that is numerical noise and can come out negative. The situation is
#' reported in a message. It arises with dichotomous items (the ladder above
#' avoids the two rungs where it is structural) and, for instance, with
#' \code{group_partial} specifications that free every parameter a rung would
#' constrain.
#'
#' Failure at a rung does not end the conversation: partial invariance
#' (freeing the offending parameter) is the usual next step, for which
#' \code{group_partial} refits the whole ladder with named parameters free,
#' and the \code{"fits"} attribute gives the fitted objects to
#' \code{\link[lavaan]{lavTestScore}} for a score test of which constraint
#' is doing the damage. This function deliberately reports the standard
#' ladder rather than automating modification searches.
#'
#' A rung whose likelihood has no maximum, flat or still rising as an
#' error variance in some group falls below zero, has no estimates to
#' compare, only the point where \pkg{lavaan}'s optimizer stopped; the
#' function checks each rung with a negative error variance as
#' \code{\link{fa_bifactor}} describes and stops with an error of class
#' \code{dmar_no_maximum} that names the rung, the group, and the item and
#' prints the \pkg{lavaan} code of that rung's fit, where
#' \code{bounds = "pos.var"} (also accepted through \code{\dots} here) or a
#' different model are the choices. Three items whose correlations in one
#' group have a negative product give such a rung when
#' \code{std.lv = TRUE} is passed through \code{\dots}; under the default
#' scaling by a marker item, \pkg{lavaan}'s optimizer does not converge on
#' such data. A negative error variance at a maximum of the likelihood is
#' a finite estimate, and its rung is reported.
#'
#' A rung that \pkg{lavaan}'s optimizer does not converge on has no
#' estimates to report, and the function stops with an error of class
#' \code{dmar_not_converged}, the class \code{\link{fa_confirmatory}}
#' gives such a fit, whose message names the rung. The condition carries
#' the rung as \code{level} and the \pkg{lavaan} code that reproduces its
#' fit as \code{syntax} (the model) and \code{lavaan_call} (the call, with
#' \code{my_data} in place of the data).
#'
#' @return A tidy wide \code{data.frame} (class \code{dmar_tbl}) with one
#'   row per fitted level and columns \code{level} (label),
#'   \code{chi_square}, \code{df}, \code{p_chi_square} (exact-fit test),
#'   \code{cfi}, \code{rmsea}, and, from the second row on, the
#'   step-comparison columns \code{delta_chi_square}, \code{delta_df},
#'   \code{p_value} (the likelihood ratio test against the previous rung),
#'   \code{delta_cfi}, and \code{delta_rmsea}. \code{delta_chi_square} and
#'   \code{p_value} are \code{NA} wherever \code{delta_df} is zero, since a
#'   difference test on zero degrees of freedom tests nothing. Attributes
#'   carry the non-numeric information: \code{"fits"}, the named list of fitted
#'   \pkg{lavaan} objects, one per level; \code{"estimator"}, the estimator
#'   actually used; \code{"ordered"}, \code{TRUE} when any indicator was
#'   declared ordered; \code{"test"}, naming the chi square difference test
#'   used (\code{NA} when a single rung was fit and nothing was compared,
#'   and for a ULS ladder of continuous items without a corrected
#'   statistic, which has no difference test);
#'   \code{"fit_indices"}, the version of the \code{cfi} and \code{rmsea}
#'   columns at every rung, \code{"standard"}, \code{"robust"}, or
#'   \code{"scaled"} (\code{"robust"} is the corrected indices of Zhang and
#'   Savalei, 2023, under FIML and two-stage estimation, and
#'   \code{"scaled"} the fallback of a ladder of ordered items, or fit by
#'   FIML or in two stages, with a scaled statistic whose robust indices
#'   \pkg{lavaan} did not compute at every rung);
#'   \code{"fit_columns"}, the lines printed under the table that name that
#'   version for a ladder fit by FIML, one that fell back from the robust
#'   indices, or one whose robust indices DMAR computed, and that say why
#'   columns of a ULS ladder of continuous items are \code{NA}, or why the
#'   \code{rmsea} of a rung with zero degrees of freedom that does not fit
#'   exactly is (\code{NULL} otherwise; a two-stage ladder names the
#'   version in its own lines); \code{"model_test"}, naming the statistic in the
#'   \code{chi_square} column at every rung (\code{"standard"},
#'   \code{"scaled"}, or under the two-stage options the two-stage residual
#'   based or scaled statistic); for a
#'   two-stage ladder, \code{"two_stage"}, a list with the method, the
#'   lines printed under the table, and, in \code{rows}, the cases the
#'   first stage read in each group and \pkg{lavaan}'s count of them when
#'   an auxiliary variable was named, how many were observed on an
#'   auxiliary variable alone, whether the rungs were computed again with
#'   all of them (\code{corrected}), and whether the estimates of a rung
#'   that constrains across groups were too (\code{reestimated}); and
#'   \code{"model"}, the \pkg{lavaan} syntax that was fitted.
#'
#' @references
#' Brosseau-Liard, P. E., & Savalei, V. (2014). Adjusting incremental fit
#'   indices for nonnormality. \emph{Multivariate Behavioral Research,
#'   49}(5), 460--470. \doi{10.1080/00273171.2014.933697}
#'
#' Brosseau-Liard, P. E., Savalei, V., & Li, L. (2012). An investigation of
#'   the sample performance of two nonnormality corrections for RMSEA.
#'   \emph{Multivariate Behavioral Research, 47}(6), 904--930.
#'   \doi{10.1080/00273171.2012.715252}
#'
#' Browne, M. W. (1974). Generalized least squares estimators in the
#'   analysis of covariance structures. \emph{South African Statistical
#'   Journal, 8}, 1--24.
#'
#' Browne, M. W. (1984). Asymptotically distribution-free methods for the
#'   analysis of covariance structures. \emph{British Journal of
#'   Mathematical and Statistical Psychology, 37}(1), 62--83.
#'
#' Byrne, B. M., Shavelson, R. J., & Muthén, B. (1989). Testing for the
#'   equivalence of factor covariance and mean structures: The issue of
#'   partial measurement invariance. \emph{Psychological Bulletin, 105}(3),
#'   456--466.
#'
#' Cheung, G. W., & Rensvold, R. B. (2002). Evaluating goodness-of-fit
#'   indexes for testing measurement invariance. \emph{Structural Equation
#'   Modeling, 9}(2), 233--255. \doi{10.1207/S15328007SEM0902_5}
#'
#' Graham, J. W. (2003). Adding missing-data-relevant variables to
#'   FIML-based structural equation models. \emph{Structural Equation
#'   Modeling, 10}(1), 80--100. \doi{10.1207/S15328007SEM1001_4}
#'
#' Meredith, W. (1993). Measurement invariance, factor analysis and
#'   factorial invariance. \emph{Psychometrika, 58}(4), 525--543.
#'   \doi{10.1007/BF02294825}
#'
#' Millsap, R. E. (2011). \emph{Statistical approaches to measurement
#'   invariance}. Routledge.
#'
#' Millsap, R. E., & Yun-Tein, J. (2004). Assessing factorial invariance in
#'   ordered-categorical measures. \emph{Multivariate Behavioral Research,
#'   39}(3), 479--515. \doi{10.1207/s15327906mbr3903_4}
#'
#' Muthén, B. (1984). A general structural equation model with dichotomous,
#'   ordered categorical, and continuous latent variable indicators.
#'   \emph{Psychometrika, 49}(1), 115--132.
#'
#' Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square test
#'   statistic for moment structure analysis. \emph{Psychometrika, 66}(4),
#'   507--514. \doi{10.1007/BF02296192}
#'
#' Savalei, V. (2018). On the computation of the RMSEA and CFI from the
#'   mean-and-variance corrected test statistic with nonnormal data in SEM.
#'   \emph{Multivariate Behavioral Research, 53}(3), 419--429.
#'   \doi{10.1080/00273171.2018.1455142}
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
#' Wu, H., & Estabrook, R. (2016). Identification of confirmatory factor
#'   analysis models of different levels of invariance for ordered
#'   categorical outcomes. \emph{Psychometrika, 81}(4), 1014--1045.
#'   \doi{10.1007/s11336-016-9506-0}
#'
#' Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods for
#'   mean and covariance structure analysis with nonnormal missing data.
#'   \emph{Sociological Methodology, 30}, 165--200.
#'
#' Zhang, X., & Savalei, V. (2023). New computations for RMSEA and CFI
#'   following FIML and TS estimation with missing data.
#'   \emph{Psychological Methods, 28}(2), 263--283.
#'   \doi{10.1037/met0000445}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cfa_1}} and \code{\link{fa_confirmatory}} for the
#'   single-group measurement models; \code{\link{reliability_omega}} for
#'   the reliability of the composite the model justifies;
#'   \code{\link{compare_cov_structures}} for covariance-structure
#'   comparisons outside the factor model;
#'   \code{\link[lavaan]{lavTestScore}} for the score test that localizes a
#'   failed rung.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords models multivariate
#'
#' @examples
#' # Do the two schools measure the verbal and the reasoning construct
#' # in the same way? (Holzinger & Swineford, bundled.) The measurement
#' # model is a named list of factors; for a single factor, name its
#' # indicators with items = instead, or give lavaan syntax, as the last
#' # example does.
#' data(holzinger_swineford)
#' hs_factors <- list(
#'   verbal    = c("t6_paragraph_comprehension", "t7_sentence",
#'                 "t9_word_meaning"),
#'   deduction = c("t20_deduction", "t22_problem_reasoning",
#'                 "t23_series_completion"))
#'
#' # The whole ladder, configural through metric, scalar, and strict.
#' # Configural invariance asks whether the same pattern of loadings
#' # holds at both schools; metric adds the constraint that the loadings
#' # are equal across schools, and it is the rung that has to hold before
#' # a relation involving one of these factors is compared across them;
#' # scalar equates the intercepts as well, which comparing factor means
#' # requires; strict equates the residual variances. Each row from the
#' # second on tests its rung against the one below it.
#' mi <- measurement_invariance(holzinger_swineford, hs_factors,
#'                              group = "school")
#' mi
#'
#' # The fitted models travel with the table, so localizing a failed rung
#' # costs no refitting.
#' names(attr(mi, "fits"))
#'
#' # The broom verbs: one row per rung of the ladder, and the model-level
#' # summary (estimator, test flavor, fit index flavor).
#' generics::tidy(mi)
#' generics::glance(mi)
#'
#' # Partial invariance frees one loading across the schools at every
#' # rung (Byrne, Shavelson, & Muthén, 1989), so the metric rung costs
#' # one degree of freedom fewer than it did above. Naming a leading
#' # subset of the rungs in levels = stops the ladder there, which keeps
#' # this and the next two examples quick; a reported analysis climbs
#' # the whole ladder.
#' measurement_invariance(holzinger_swineford, hs_factors,
#'                        group = "school",
#'                        group_partial = "verbal =~ t7_sentence",
#'                        levels = c("configural", "metric"))
#'
#' # With the indicators coded as ordered categories, ordered = TRUE puts
#' # the thresholds rung first, because thresholds rather than intercepts
#' # carry the location information there, and makes the rung-to-rung
#' # tests the scaled difference tests; the "test" attribute names the
#' # test that was used. The ladder stops at the thresholds rung here;
#' # metric, scalar, and strict follow it. lavaan computes no robust fit
#' # indices for the thresholds rung, so both rows report the scaled ones,
#' # as the line under the table says.
#' hs_ordered <- holzinger_swineford
#' for (item in unlist(hs_factors, use.names = FALSE)) {
#'   hs_ordered[[item]] <- as.integer(cut(
#'     holzinger_swineford[[item]],
#'     breaks = quantile(holzinger_swineford[[item]],
#'                       c(0, .25, .5, .75, 1)),
#'     include.lowest = TRUE))
#' }
#' mi_ordered <- measurement_invariance(hs_ordered, hs_factors,
#'                                      group = "school", ordered = TRUE,
#'                                      levels = c("configural",
#'                                                 "thresholds"))
#' mi_ordered
#' attr(mi_ordered, "test")
#'
#' # The measurement model can also be given as lavaan syntax, and
#' # missing = "fiml", the default, fits the ladder by full information
#' # maximum likelihood when continuous items have incomplete cases; here
#' # twenty scores on one item are set to missing first. Naming it states
#' # the choice in the call. With three indicators
#' # the single factor is just identified in each group, so the
#' # configural row is saturated (zero degrees of freedom, exact fit)
#' # and the metric row carries the test. Under FIML the cfi and rmsea
#' # columns are the corrected indices of Zhang and Savalei (2023), as the
#' # line under the table says.
#' hs_missing <- holzinger_swineford
#' set.seed(113)
#' hs_missing$t7_sentence[sample(nrow(hs_missing), 20)] <- NA
#' measurement_invariance(
#'   hs_missing,
#'   model = "verbal =~ t6_paragraph_comprehension + t7_sentence +
#'            t9_word_meaning",
#'   group = "school", missing = "fiml",
#'   levels = c("configural", "metric"))
#'
#' @export
measurement_invariance <- function(data, model = NULL, group,
                                   items = NULL,
                                   levels = NULL,
                                   ordered = NULL,
                                   estimator = "ML",
                                   missing = "fiml",
                                   group_partial = NULL,
                                   parameterization = c("delta", "theta"),
                                   ...) {
  # Whether 'missing' was left at its default, read before anything else:
  # the default resolves by rule (.measurement_invariance_resolve_missing()).
  missing_default <- base::missing(missing)
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    .dmar_refuse(c(
      "Package 'lavaan' is required for measurement_invariance(). ",
      "Install it with install.packages(\"lavaan\")."),
      "package_required", package = "lavaan")
  }
  if (!is.data.frame(data)) {
    .dmar_refuse_argument("'data' must be a data.frame.", arg = "data")
  }
  parameterization <- match.arg(parameterization)

  # ---- the measurement model: syntax, a named list, or one factor ------
  if (is.null(model) == is.null(items)) {
    .dmar_refuse_argument(c(
      "Supply exactly one of 'model' (lavaan syntax or a named list of ",
      "factors) and 'items' (the indicators of a single factor)."),
      arg = c("model", "items"))
  }
  # Items that cannot be read from 'data' are refused with the class the
  # factor analysis functions give the mistake, dmar_items_invalid, and
  # the reasons of .fa_items_check() (R/factor_analysis_internals.R), so
  # that a program reads a misnamed item from the class rather than from
  # the wording. The checks run in that function's order (form, unknown,
  # duplicate, count), and a repeated item is refused with the message
  # measurement_alignment() gives it.
  if (!is.null(items)) {
    failed <- .fa_items_check(items, names(data))
    if (!is.null(failed)) {
      .fa_items_refuse(if (identical(failed$reason, "duplicate")) {
        "'items' must not repeat an item."
      } else {
        "'items' must name three or more columns of 'data'."
      }, failed$reason, failed$items)
    }
    syntax <- paste("f =~", paste(items, collapse = " + "))
  } else if (is.list(model)) {
    syntax <- .measurement_invariance_syntax(model, data)
  } else if (is.character(model) && length(model) >= 1L &&
             !anyNA(model)) {
    syntax <- paste(model, collapse = "\n")
    if (!grepl("=~", syntax, fixed = TRUE)) {
      .dmar_refuse_argument(
        "'model' must contain at least one measurement relation (=~).",
        arg = "model")
    }
  } else {
    .dmar_refuse_argument(c(
      "'model' must be lavaan model syntax (a character string) or a ",
      "named list: list(factor_name = c(\"item_1\", \"item_2\", ...))."),
      arg = "model")
  }

  # Indicators of the fitted model, needed for ordered = TRUE and to check
  # that every indicator is actually in the data. One the data lack is
  # refused as a misnamed item (.fa_items_absent(), reason "unknown").
  parsed <- lavaan::lavaanify(syntax)
  latent <- unique(parsed$lhs[parsed$op == "=~"])
  indicators <- setdiff(unique(parsed$rhs[parsed$op == "=~"]), latent)
  absent <- setdiff(indicators, names(data))
  if (length(absent)) {
    .fa_items_absent(c("Indicators not found in 'data': ",
                       paste(absent, collapse = ", "), "."), absent)
  }

  # ---- the grouping variable -------------------------------------------
  if (!is.character(group) || length(group) != 1L ||
      !(group %in% names(data))) {
    .dmar_refuse_argument(
      "'group' must name one column of 'data'.", arg = c("group", "data"))
  }
  if (length(unique(stats::na.omit(data[[group]]))) < 2L) {
    .dmar_refuse_argument(
      "'group' must have at least two groups.", arg = "group")
  }

  # ---- estimator, missing data, ordered indicators ---------------------
  estimator_choices <- c("ML", "MLM", "MLMV", "MLR", "GLS", "WLS", "DWLS",
                         "WLSM", "WLSMV", "ULS", "ULSM", "ULSMV")
  missing_choices <- c("fiml", "ml", "listwise", "pairwise",
                       .TWO_STAGE_VALUES)
  if (!is.character(estimator) || length(estimator) != 1L ||
      !(estimator %in% estimator_choices)) {
    .dmar_refuse_argument(c(
      "'estimator' must be one of ",
      paste(dQuote(estimator_choices), collapse = ", "), "."),
      arg = "estimator")
  }
  if (!is.character(missing) || length(missing) != 1L ||
      !(missing %in% missing_choices)) {
    .dmar_refuse_argument(c(
      "'missing' must be one of ",
      paste(dQuote(missing_choices), collapse = ", "), "."), arg = "missing")
  }
  missing <- .measurement_invariance_resolve_missing(
    missing, missing_default, data, lavaan::lavNames(parsed, "ov"), group,
    ordered = isTRUE(ordered) ||
      (is.character(ordered) && length(ordered) > 0L),
    estimator = estimator, se = list(...)[["se"]],
    test = list(...)[["test"]])
  # lavaan's value for a two-stage request, NULL otherwise.
  ts_missing <- .ts_lavaan_missing(missing)
  two_stage <- !is.null(ts_missing)

  ordered_items <- NULL
  if (!is.null(ordered)) {
    if (isTRUE(ordered)) {
      ordered_items <- indicators
    } else if (is.character(ordered) && length(ordered)) {
      unknown <- setdiff(ordered, indicators)
      if (length(unknown)) {
        .dmar_refuse_argument(c(
          "'ordered' names variables that are not indicators of the ",
          "model: ", paste(unknown, collapse = ", "), "."), arg = "ordered")
      }
      ordered_items <- unique(ordered)
    } else {
      .dmar_refuse_argument(c(
        "'ordered' must be TRUE (all indicators) or a character vector ",
        "of item names."), arg = "ordered")
    }
    if (missing %in% c("ml", "fiml")) {
      .dmar_refuse_argument(c(
        "Full information maximum likelihood is not available with the ",
        "categorical estimator; with ordered items use ",
        "missing = \"pairwise\" or missing = \"listwise\"."), arg = NULL)
    }
    if (two_stage) {
      .dmar_refuse_argument(c(
        "Two-stage estimation is not available with ordered items: its ",
        "first stage estimates the means and covariances of continuous ",
        "variables under the normal model (lavaan refuses it in the ",
        "categorical setting). Use missing = \"pairwise\" or ",
        "missing = \"listwise\"."), arg = NULL)
    }
    if (estimator %in% c("ML", "MLM", "MLMV", "MLR", "GLS")) {
      message("Ordered items declared: switching to estimator = \"WLSMV\".")
      estimator <- "WLSMV"
    } else if (estimator %in% c("DWLS", "ULS")) {
      # Neither carries a robust test, so the rung-to-rung difference test
      # would be the unscaled one, which is not asymptotically valid for
      # ordered data. The promotion changes only the test statistic: in
      # lavaan "WLSMV" is "DWLS" and "ULSMV" is "ULS" with the mean and
      # variance adjustment, so the discrepancy function, and with it every
      # parameter estimate, is untouched. "WLS" is deliberately left alone
      # (see the estimator documentation).
      promoted <- c(DWLS = "WLSMV", ULS = "ULSMV")[[estimator]]
      message("Ordered items declared: switching to estimator = \"",
              promoted, "\". That is the same ", estimator,
              " fit with the mean and variance adjusted test, which the ",
              "chi square difference test between rungs needs to be ",
              "asymptotically valid with ordered data; the parameter ",
              "estimates are unchanged.")
      estimator <- promoted
    }
  }
  is_ordered <- !is.null(ordered_items)

  # Two-stage estimation fits the EM moments by maximum likelihood, and its
  # standard errors are the two-stage sandwich (robust to nonnormality
  # under "robust_two_stage"); lavaan would override another estimator or
  # a robust 'se' with a warning, which here is an error naming the option
  # that does what was asked, as in fa_confirmatory().
  dots <- list(...)
  if (two_stage) {
    if (estimator != "ML") {
      .dmar_refuse_argument(c(
        "Two-stage estimation fits the EM estimates of the means and ",
        "covariances by maximum likelihood (estimator = \"ML\"). ",
        if (estimator == "MLR") {
          paste0("For standard errors and tests robust to nonnormality ",
                 "use missing = \"robust_two_stage\" rather than ",
                 "estimator = \"MLR\".")
        } else {
          paste0("estimator = \"", estimator, "\" is not available ",
                 "with it.")
        }), arg = NULL)
    }
    if (!is.null(dots$se)) {
      if (!dots$se %in% c("standard", "none", ts_missing)) {
        .dmar_refuse_argument(c(
          "Under two-stage estimation the standard errors are the ",
          "two-stage ones, so 'se' must be \"standard\" (or \"none\").",
          if (dots$se %in% c("robust.sem", "robust.huber.white")) {
            paste0(" For standard errors robust to nonnormality use ",
                   "missing = \"robust_two_stage\".")
          }), arg = NULL)
      }
      dots$se <- ts_missing
    }
  }
  # The auxiliary variables passed through ... as lavaan's aux. Under FIML
  # lavaan adds each one to the model as a saturated correlate (Graham,
  # 2003), with a mean, a variance, and covariances with every observed
  # variable as ordinary parameters of each group, so group.equal reaches
  # them: the intercepts rung would equate the auxiliary's means across
  # groups and the residuals rung its variances, and the rung tests would
  # then test whether the auxiliary is distributed alike in the groups,
  # which says nothing about the measurement model. They are freed at every
  # rung through group.partial, beside any terms the caller frees. (No rung
  # equates the covariances, since "residual.covariances" is not on the
  # ladder.) Through version 1.0.1 they were equated: on six Holzinger and
  # Swineford tests by school, with 60 sentence completion scores missing
  # and general information as the auxiliary, whose means are 44.8 and 36.6
  # in the two schools, the scalar rung gave a difference of 57.59 on 5 df,
  # against 24.03 on 4 without the auxiliary. Under two-stage estimation
  # the auxiliaries enter the first stage only and carry no parameter.
  aux_dots <- .ts_dots_aux(dots)
  # Under FIML, by default or by request: an explicit request that lavaan
  # would fit by listwise deletion instead is refused, and so are data in
  # which a group never observes a variable of the model (or an auxiliary
  # variable), or never observes two together, which lavaan cannot fit.
  .measurement_invariance_check_fiml_request(missing, missing_default,
                                             estimator, dots[["se"]],
                                             dots[["test"]])
  if (missing %in% c("ml", "fiml")) {
    .measurement_invariance_check_coverage(
      data, lavaan::lavNames(parsed, "ov"), aux_dots, group)
  }
  aux_partial <- if (missing %in% c("ml", "fiml") && length(aux_dots)) {
    c(paste(aux_dots, "~ 1"), paste(aux_dots, "~~", aux_dots))
  }
  # Under two-stage estimation with auxiliary variables lavaan warns, at
  # every rung, that a case observed on an auxiliary variable alone is
  # empty and "will be ignored". Its first stage reads such a case, and
  # every rung is computed again with it counted below (.ts_all_rows()),
  # so the warning is wrong about it; it is muffled when no case is empty
  # of the auxiliary variables as well, and .ts_aux_only_note() says what
  # is done instead.
  ts_quiet_empty <- two_stage && length(aux_dots) > 0L &&
    !any(rowSums(!is.na(data[, intersect(
      c(lavaan::lavNames(parsed, "ov"), aux_dots), names(data)),
      drop = FALSE])) == 0L)

  # A dichotomous indicator carries one threshold, and that threshold is not
  # separately identified from the intercept of its underlying response, nor
  # is its residual variance a free parameter (Millsap & Yun-Tein, 2004; Wu &
  # Estabrook, 2016). Counting the observed categories here is what keeps the
  # ladder from including rungs that cannot be tested. The strict rung goes
  # only when no indicator is left to carry a free residual variance, since a
  # continuous indicator alongside dichotomous ones still has one.
  binary_items <- character(0)
  all_binary <- FALSE
  drop_strict <- FALSE
  if (is_ordered) {
    n_categories <- .measurement_invariance_n_categories(data, ordered_items)
    binary_items <- ordered_items[n_categories <= 2L]
    all_binary <- length(binary_items) == length(ordered_items)
    drop_strict <- all_binary && setequal(ordered_items, indicators)
  }

  # An empty set of terms, in either spelling below (NULL, character(0), or
  # "", lavaan's own default), frees nothing: lavaan reads each as no
  # exemption from group.equal, so every rung imposes its full set of
  # constraints, the full invariance ladder (Meredith, 1993), and partial
  # invariance (Byrne, Shavelson, & Muthén, 1989) begins only when a term
  # is named. Through 1.0.1 an empty group.partial in ... reached lavaan
  # and was read so, and the formal argument refused character(0).
  if (!is.null(group_partial) && !length(group_partial)) group_partial <- NULL
  if (!is.null(group_partial) &&
      (!is.character(group_partial) || anyNA(group_partial))) {
    .dmar_refuse_argument(c(
      "'group_partial' must be a character vector of parameter labels ",
      "or lavaan parameter specifications, such as ",
      "\"visual =~ item_2\"."), arg = "group_partial")
  }
  # lavaan's own spelling, group.partial passed through ..., frees its terms
  # at every rung as group_partial does. The two are merged, with the terms
  # that free an auxiliary variable under FIML (aux_partial above), into
  # the one group.partial each rung passes to lavaan, which stops on a
  # second with an error that a name is duplicated (through version 1.0.1
  # a call that gave both spellings stopped so). group.equal is what the
  # ladder sets rung by rung, so it is refused in ...: through version
  # 1.0.1 it constrained the configural rung, which was then reported as
  # configural, and stopped every later rung with the same error.
  dot_names <- names(dots)
  if (any(dot_names %in% c("group.equal", "group_equal"))) {
    .dmar_refuse_argument(c(
      "'group.equal' is set by the ladder at every rung, so it cannot ",
      "be passed through '...'. Choose the rungs with 'levels', and ",
      "leave parameters free across groups with 'group_partial'."),
      arg = c("levels", "group_partial"))
  }
  in_dots <- dot_names %in% "group.partial"
  if (any(in_dots)) {
    dots_partial <- unlist(dots[in_dots], use.names = FALSE)
    dots <- dots[!in_dots]
    if (length(dots_partial)) {
      if (!is.character(dots_partial) || anyNA(dots_partial)) {
        .dmar_refuse_argument(c(
          "'group.partial' must be a character vector of parameter ",
          "labels or lavaan parameter specifications, such as ",
          "\"visual =~ item_2\"."), arg = NULL)
      }
      group_partial <- unique(c(group_partial, dots_partial))
    }
  }
  group_partial <- group_partial[nzchar(group_partial)]
  if (!length(group_partial)) group_partial <- NULL

  # ---- the ladder in force ---------------------------------------------
  # With every declared item dichotomous the thresholds rung is folded into
  # metric (which constrains thresholds and loadings together), and the
  # strict rung goes too when no indicator is left to carry a free residual
  # variance. Each dropped rung would have cost zero degrees of freedom.
  ladder <- if (!is_ordered) {
    c("configural", "metric", "scalar", "strict")
  } else {
    rungs <- c("configural", "thresholds", "metric", "scalar", "strict")
    if (all_binary) rungs <- setdiff(rungs, "thresholds")
    if (drop_strict) rungs <- setdiff(rungs, "strict")
    rungs
  }
  ladder_label <- if (!is_ordered) {
    "continuous"
  } else if (all_binary) {
    "dichotomous ordered"
  } else {
    "ordered"
  }
  if (all_binary) {
    message("Every declared ordered indicator is dichotomous. A threshold ",
            "and the intercept of the underlying response are then not ",
            "separately identified (Millsap & Yun-Tein, 2004; Wu & ",
            "Estabrook, 2016), so the thresholds rung is folded into metric, ",
            "which constrains thresholds and loadings together. ",
            if (drop_strict) {
              paste0("The strict rung is dropped as well, since a ",
                     "dichotomous indicator has no free residual variance to ",
                     "constrain. ")
            } else {
              paste0("The strict rung is kept, since it still constrains the ",
                     "residual variances of the indicators that were not ",
                     "declared ordered. ")
            },
            "The ladder is ", paste(ladder, collapse = ", "), ".")
  } else if (length(binary_items)) {
    message(length(binary_items), " of ", length(ordered_items),
            " declared ordered indicators are dichotomous (",
            paste(binary_items, collapse = ", "),
            "). The full ordered ladder is kept, because the polytomous ",
            "items carry testable thresholds and residual variances, but ",
            "the thresholds and strict rungs test nothing about the ",
            "dichotomous ones.")
  }
  constraints <- if (is_ordered) {
    list(configural = character(0),
         thresholds = "thresholds",
         metric     = c("thresholds", "loadings"),
         scalar     = c("thresholds", "loadings", "intercepts"),
         strict     = c("thresholds", "loadings", "intercepts",
                        "residuals"))
  } else {
    list(configural = character(0),
         metric     = "loadings",
         scalar     = c("loadings", "intercepts"),
         strict     = c("loadings", "intercepts", "residuals"))
  }
  if (is.null(levels)) {
    levels <- ladder
  } else {
    if (!is.character(levels) || !length(levels) || anyNA(levels)) {
      .dmar_refuse_argument(
        "'levels' must be a character vector of rung names.", arg = "levels")
    }
    unknown <- setdiff(levels, ladder)
    if (length(unknown)) {
      .dmar_refuse_argument(c(
        "'levels' names rungs that are not on the ", ladder_label,
        " ladder: ", paste(unknown, collapse = ", "),
        ". The ladder in force is ", paste(ladder, collapse = ", "), ".",
        if (all_binary)
          paste0(" With every declared ordered indicator dichotomous, a",
                 " rung missing from that list would impose no",
                 " restriction; see ?measurement_invariance.")), arg = "levels")
    }
    levels <- ladder[ladder %in% levels]
    if (!identical(levels, ladder[seq_along(levels)])) {
      .dmar_refuse_argument(c(
        "'levels' must be a leading subset of the ladder in force: ",
        paste(ladder, collapse = ", "), "."), arg = "levels")
    }
  }

  # Residual variances are estimated only under the theta parameterization,
  # so a strict rung on ordered items would otherwise cost zero degrees of
  # freedom and test nothing.
  if (is_ordered && "strict" %in% levels && parameterization == "delta") {
    message("A strict rung with ordered items requires the theta ",
            "parameterization (residual variances are not free parameters ",
            "under delta): switching to parameterization = \"theta\".")
    parameterization <- "theta"
  }

  # ---- fit the ladder ---------------------------------------------------
  rung_args <- function(lv) {
    args <- list(model = syntax, data = data, group = group,
                 estimator = estimator,
                 missing = if (two_stage) ts_missing else missing,
                 meanstructure = TRUE)
    if (is_ordered) {
      args$ordered <- ordered_items
      args$parameterization <- parameterization
    }
    if (length(constraints[[lv]])) args$group.equal <- constraints[[lv]]
    partial <- c(group_partial, aux_partial)
    if (length(partial)) args$group.partial <- partial
    c(args, dots)
  }
  fits <- lapply(levels, function(lv) {
    withCallingHandlers(
      do.call(lavaan::cfa, rung_args(lv)),
      warning = function(cond) {
        if (ts_quiet_empty &&
            grepl("some cases are empty and will be ignored",
                  conditionMessage(cond), fixed = TRUE))
          invokeRestart("muffleWarning")
      })
  })
  names(fits) <- levels
  .measurement_invariance_check_fiml_fits(fits, missing)
  # ---- a likelihood with no maximum ------------------------------------
  # Each rung, before anything is read from it: when an error variance is
  # negative in some group, whether the likelihood has a maximum at all
  # (.lavaan_no_maximum(), R/no_maximum_internals.R). A rung on a ridge,
  # where that variance runs off below zero with the likelihood flat or
  # still rising and lavaan's optimizer stopped where its steps stopped
  # gaining, is refused, with the rung and the group named and the lavaan
  # code that reproduces it; a finite negative estimate is not. The check
  # does not apply to the WLSMV ladder of ordered items, or to least
  # squares, which have no likelihood. A rung that did not converge is
  # refused first, rung by rung, with the class fa_confirmatory() gives
  # such a fit (.measurement_invariance_not_converged()); through 1.0.x and
  # the first versions of 1.1.0 lavaan's own error from its fit measures
  # ("fit measures not available if model did not converge") ended the
  # call, without naming the rung.
  for (lv in levels) {
    if (!isTRUE(lavaan::lavInspect(fits[[lv]], "converged"))) {
      .measurement_invariance_not_converged(
        lv, syntax, .no_max_call_text("lavaan::cfa", rung_args(lv)))
    }
    no_max <- .lavaan_no_maximum(fits[[lv]])
    if (!is.null(no_max)) {
      .no_max_refuse("measurement_invariance", no_max, syntax = syntax,
                     lavaan_call = .no_max_call_text("lavaan::cfa",
                                                     rung_args(lv)),
                     where = paste("the", lv, "rung"),
                     # The model may be the user's own syntax, whose ":="
                     # lines define the user's quantities, so it is printed
                     # whole; the syntax DMAR builds from items has none.
                     defined = "keep")
    }
  }
  # A case observed on an auxiliary variable but on none of the model's
  # variables is read by the first stage and left out of lavaan's N, so
  # every rung is computed again with its N set to the cases the first
  # stage read (.ts_all_rows(), which fa_confirmatory() applies too): the statistics,
  # the indices, the difference tests, and the fits returned in the "fits"
  # attribute use every case. The configural rung keeps lavaan's
  # estimates, which do not depend on the counts when no parameter is
  # constrained across groups; the rungs that constrain across groups are
  # estimated again, since the discrepancy function weights each group by
  # its count and the groups' counts can grow in different ratios. Through
  # the first versions of 1.1.0 such cases were set aside, which
  # conditioned the estimates on an observed item.
  ts_rows <- NULL
  if (two_stage) {
    all_rows <- lapply(fits, .ts_all_rows)
    # The rungs read the same data, so they have the same cases. The
    # rungs are compared with one another, so they are all computed again
    # or, if one could not be, none is.
    ts_rows <- all_rows[[1L]][c("read", "nobs", "extra", "corrected",
                                "reestimated")]
    if (ts_rows$extra > 0L) {
      ts_rows$corrected <- all(vapply(all_rows, function(r)
        isTRUE(r$corrected), logical(1L)))
      ts_rows$reestimated <- any(vapply(all_rows, function(r)
        isTRUE(r$reestimated), logical(1L)))
      if (ts_rows$corrected) {
        fits <- lapply(all_rows, `[[`, "fit")
        names(fits) <- levels
        .ts_aux_only_note(
          "measurement_invariance", ts_rows$extra, sum(ts_rows$read),
          sum(ts_rows$nobs),
          what = paste0(
            if (ts_rows$reestimated) {
              paste0("the estimates of the rungs that constrain across ",
                     "groups, ")
            },
            "the test statistics, fit indices, and difference tests of ",
            "every rung are computed on all ", sum(ts_rows$read),
            " cases it reads"))
      }
      .ts_all_rows_warning("measurement_invariance", ts_rows, ts_missing)
    }
  }

  # ---- fit indices ------------------------------------------------------
  # With a robust or categorical estimator the standard statistics are not
  # the ones to report; lavaan's scaled and robust versions are. Whether the
  # fit carries such a statistic is read from its test option, setting aside
  # the statistics that are not corrections of the standard one (see
  # .lavaan_corrected_tests()).
  scaled <- length(.lavaan_corrected_tests(
    lavaan::lavInspect(fits[[1L]], "options")$test)) > 0L
  # Under FIML the cfi and rmsea columns are lavaan's robust versions, the
  # corrected indices of Zhang and Savalei (2023), as in fa_confirmatory(): they
  # estimate the values complete data would give, the values a cutoff for
  # a change in fit (Cheung & Rensvold, 2002) is a statement about, where
  # the standard indices computed from the FIML statistic change with the
  # amount of missing data. Through 1.0.1 an ML ladder under FIML
  # reported the standard ones. Whether the ladder used FIML is read from
  # the fit rather than from the argument: lavaan fits the complete cases
  # instead, with a warning, when missing = "ml" meets a test it cannot
  # combine with FIML (test = "satorra.bentler" or "mean.var.adjusted"
  # passed through ...), and the indices of such a fit are not the FIML
  # corrected ones.
  fiml <- isTRUE(lavaan::lavInspect(fits[[1L]], "options")$missing %in%
                   c("ml", "ml.x"))
  # How the cfi and rmsea columns fall back when lavaan has no robust
  # versions at every rung (.measurement_invariance_indices()): by the
  # missing data method under FIML and two-stage estimation, and otherwise
  # by the estimator, "ml" for a maximum likelihood fit of continuous items
  # (complete data or listwise deletion) and "other" for ordered items or
  # least squares.
  method <- if (fiml) {
    "fiml"
  } else if (two_stage) {
    "two_stage"
  } else if (!is_ordered && identical(
    lavaan::lavInspect(fits[[1L]], "options")$estimator, "ML")) {
    "ml"
  } else "other"
  # The measures every version could need, read once for every rung: the
  # versions are chosen from them and the columns are taken from them. One
  # version of the test and one of the indices hold every rung, so that the
  # columns of a rung and the differences between rungs compare like with
  # like (.measurement_invariance_test_columns() and
  # .measurement_invariance_indices()). Through 1.0.1 a scaled ladder read
  # the names from its first rung alone and reported the standard indices
  # where lavaan had not computed robust ones (it computes none beside test
  # = "mean.var.adjusted") while the fit_indices attribute said "robust",
  # and the standard statistic where it had no scaled one (beside test =
  # "bollen.stine") while the model_test attribute said "scaled".
  h1_issues <- lapply(fits, .fiml_fit_h1_diagnostic)
  if (any(!vapply(h1_issues, is.null, logical(1L)))) {
    .dmar_refuse_argument(c(
      "The invariance table requires unrestricted-H1 fit indices that ",
      "are unavailable for at least one rung. Inspect the constrained ",
      "models separately. An ordinary nested ML likelihood ratio can ",
      "still compare valid H0 fits on the same cases; the absolute fit ",
      "indices and their changes cannot be reported. ",
      paste(unique(unlist(h1_issues)), collapse = " ")), arg = NULL)
  }
  fms <- lapply(fits, lavaan::fitMeasures,
                .measurement_invariance_measures(scaled, method))
  # The robust indices of a maximum likelihood ladder with the older mean
  # and variance adjusted statistic (test = "mean.var.adjusted" through
  # ...), which lavaan does not compute and DMAR does (.sem_robust_mva()):
  # through 1.0.1 such a ladder reported the standard indices, labeled
  # "robust".
  mva <- lapply(fits, .sem_robust_mva)
  fms <- Map(function(fm, extra) c(fm, extra), fms, mva)
  uls_ladder <- !is_ordered && .sem_uls_continuous(fits[[1L]])
  # The rungs with no degrees of freedom and a statistic of zero
  # (.sem_exact_fit()), which fit exactly in every version
  # (see below); the line naming robust indices DMAR computed is for the
  # rungs that have some.
  saturated <- vapply(fms, .sem_exact_fit, logical(1L))
  tc <- .measurement_invariance_test_columns(fms, scaled)
  fi <- .measurement_invariance_indices(
    fms, scaled, method, dmar_robust = any(lengths(mva) > 0L & !saturated),
    allow_robust = !uls_ladder)
  wanted <- c(tc$measures, fi$measures)
  fit_indices <- fi$version
  # The line naming the version of the cfi and rmsea columns: always under
  # FIML, in place of the two-stage line when a two-stage ladder falls
  # back, and otherwise when a scaled ladder does.
  fit_columns <- if (method != "two_stage") fi$lines

  M <- t(vapply(fms, function(fm) as.numeric(fm[wanted]),
                numeric(length(wanted))))
  # A saturated rung (df = 0) reproduces the moments exactly, so its
  # statistic is zero in every version, and with it the model term of the
  # CFI, T - c df, whatever the scaling constant c (0/0 there, which is how
  # lavaan's robust CFI comes out NA or 1 by rounding): its CFI is 1 in
  # every version (and by Bentler's, 1990, convention when the baseline
  # term is zero too). Its RMSEA is 0/0, sqrt((T - c df)/(N df)), which no
  # version defines; lavaan reports 0 there in every version where it
  # computes one, its convention for a model that fits exactly, and the
  # rung keeps it, so that delta_cfi and delta_rmsea at the next rung
  # compare with exact fit. Through 1.0.1 a robust ladder whose configural
  # rung was saturated reported its CFI as NA, and the metric rung's
  # delta_cfi with it.
  M[saturated, 4L] <- 1
  M[saturated, 5L] <- 0
  # A rung with no degrees of freedom whose statistic is not zero
  # (.sem_zero_df_misfit(), as when an inequality bound binds) does not fit
  # exactly: its cfi is formed from its statistic, and its RMSEA, which
  # divides by the degrees of freedom, is undefined, not the 0 of exact fit
  # lavaan reports whenever df = 0, so it is NA, and delta_rmsea at the
  # next rung with it, as sem_fit_measures() reports the RMSEA of such a
  # fit. Through 1.0.1 (and the first version of the rule above, which
  # stopped taking such a rung to fit exactly) the column held lavaan's 0.
  zero_df_misfit <- vapply(fms, .sem_zero_df_misfit, logical(1L))
  M[zero_df_misfit, 5L] <- NA_real_
  if (any(zero_df_misfit)) {
    fit_columns <- c(fit_columns, fit_columns_zero_df =
                       .measurement_invariance_zero_df_line(
                         levels[zero_df_misfit]))
  }
  # Continuous ULS cannot use ML-derived robust indices. Scaled indices
  # remain qualified diagnostics only where their estimates and statistic
  # are invariant to a common change of units. Mean-restricted rungs and
  # uncorrected ULS discrepancies need the additional refusals below.
  uls_restricted <- character(0)
  if (uls_ladder) {
    browne <- any(startsWith(
      lavaan::lavInspect(fits[[1L]], "options")$test, "browne.residual"))
    if (fi$version != "standard" || !browne) {
      restricted <- !saturated & vapply(fits, function(x) {
        !isFALSE(.sem_means_restricted(x))
      }, logical(1L))
      no_rmsea <- !saturated & (fi$version != "scaled" | restricted)
      no_cfi <- !saturated & (fi$version != "scaled" | restricted)
      M[no_rmsea, 5L] <- NA_real_
      M[no_cfi, 4L] <- NA_real_
      uls_restricted <- levels[restricted]
      line <- .measurement_invariance_uls_line(fi$version, uls_restricted,
                                               any(no_rmsea))
    } else {
      restricted <- !saturated & vapply(fits, function(x) {
        !isFALSE(.sem_uls_estimates_units(x))
      }, logical(1L))
      M[restricted, 4:5] <- NA_real_
      uls_restricted <- levels[restricted]
      line <- .measurement_invariance_uls_line("browne", uls_restricted,
                                               any(restricted))
    }
    policy <- paste0("Continuous ULS: ML-derived robust indices are withheld ",
      "(Savalei, 2018). Scaled indices do not estimate the same population ",
      "discrepancy as ML indices, and ML cutoffs do not automatically apply.")
    fit_columns <- c(fit_columns, fit_columns_uls = paste(policy, line))
  }
  # After normal theory two-stage estimation the test of each rung is the
  # residual based statistic of Savalei and Bentler (2009), the one fa_confirmatory()
  # reports for the same model, rather than lavaan's scaled approximation.
  # The indices stay lavaan's robust ones (the corrected indices of Zhang &
  # Savalei, 2023; see .ts_indices_label()).
  if (identical(ts_missing, "two.stage")) {
    for (i in seq_along(fits)) {
      rt <- .ts_residual_test(fits[[i]])
      M[i, 1:3] <- c(rt$statistic, rt$df, rt$p_value)
    }
  }

  # ---- the rung-to-rung comparisons -------------------------------------
  k <- length(levels)
  d_chi <- d_df <- p_lrt <- d_cfi <- d_rmsea <- rep(NA_real_, k)
  # With a single rung there is nothing to compare, so no difference test
  # was performed and the attribute says so rather than naming one.
  test_label <- NA_character_
  untestable <- character(0)
  if (k > 1L) {
    for (i in 2:k) {
      # A rung that adds no degree of freedom makes lavTestLRT warn that two
      # models share their degrees of freedom. The package says what that
      # means in its own words below, so the warning is muffled here rather
      # than reaching the user twice, once cryptically. Every other lavaan
      # warning passes through untouched.
      lrt <- withCallingHandlers(
        lavaan::lavTestLRT(fits[[i - 1L]], fits[[i]]),
        warning = function(cond) {
          if (grepl("same degrees of freedom", conditionMessage(cond),
                    fixed = TRUE))
            invokeRestart("muffleWarning")
        })
      d_df[i] <- lrt[2L, "Df diff"]
      # Zero degrees of freedom means the rung imposed no restriction: the
      # two models are reparameterizations of one another, and their chi
      # square difference is numerical noise that can come out negative. A
      # difference test is reported only where one was performed.
      if (!is.na(d_df[i]) && d_df[i] <= 0) {
        untestable <- c(untestable, levels[i])
      } else {
        d_chi[i] <- lrt[2L, "Chisq diff"]
        p_lrt[i] <- lrt[2L, "Pr(>Chisq)"]
      }
      d_cfi[i]   <- M[i, 4L] - M[i - 1L, 4L]
      d_rmsea[i] <- M[i, 5L] - M[i - 1L, 5L]
      # lavTestLRT() scales the difference only when the rungs carry a
      # scaled statistic, which a fit with a bootstrap test (test =
      # "bollen.stine") does not, though its test option lists one.
      if (i == 2L) {
        test_label <- .measurement_invariance_test_label(
          lrt, tc$version == "scaled")
      }
    }
    # lavaan's difference test of a ULS ladder of continuous items without
    # a corrected statistic refers the difference of the ULS statistics,
    # which are in the units of the squared covariances and have no chi
    # square distribution, to the chi square: through 1.0.1 a metric rung
    # of the first example under estimator = "ULS" had a delta_chi_square
    # of 91712 (p < .0001), and of 9.17 (p = .057) with every score divided
    # by 10, while the attribute named a scaled difference test. There is
    # no difference test to report, and the line under the table says so.
    if (uls_ladder && tc$version == "standard") {
      d_chi[] <- NA_real_
      p_lrt[] <- NA_real_
      test_label <- NA_character_
      fit_columns <- c(fit_columns, fit_columns_uls_test = paste0(
        "delta_chi_square, p_value: lavaan's difference test of a ULS ",
        "ladder of continuous items without a corrected statistic refers ",
        "the difference of the ULS statistics, which are in the units of ",
        "the squared covariances and are not chi square distributed, to the ",
        "chi square distribution, so there is none; estimator = \"ULSM\" or ",
        "\"ULSMV\" gives a scaled difference test, and estimator = \"MLR\" ",
        "one that does not depend on the units"))
    }
    if (length(untestable)) {
      message("The constraints at the ", paste(untestable, collapse = ", "),
              if (length(untestable) > 1L) " rungs" else " rung",
              " cost no degrees of freedom, so nothing was tested there: ",
              "delta_chi_square and p_value are NA rather than a difference ",
              "test on zero degrees of freedom.")
    }
  }

  out <- data.frame(
    level            = levels,
    chi_square       = M[, 1L],
    df               = M[, 2L],
    p_chi_square     = M[, 3L],
    cfi              = M[, 4L],
    rmsea            = M[, 5L],
    delta_chi_square = d_chi,
    delta_df         = d_df,
    p_value          = p_lrt,
    delta_cfi        = d_cfi,
    delta_rmsea      = d_rmsea,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  res <- .as_dmar_tbl(out, subclass = "dmar_measurement_invariance")
  attr(res, "fits") <- fits
  attr(res, "estimator") <- estimator
  attr(res, "ordered") <- is_ordered
  attr(res, "test") <- test_label
  attr(res, "fit_indices") <- fit_indices
  attr(res, "fit_columns") <- fit_columns
  attr(res, "model_test") <- if (two_stage) {
    .ts_test_label(ts_missing)
  } else {
    tc$version
  }
  if (two_stage) {
    ts_footer <- .ts_footer(ts_missing, "invariance", rows = ts_rows)
    if (length(fi$lines)) ts_footer[["two_stage_indices"]] <- fi$lines[[1L]]
    attr(res, "two_stage") <- list(missing = missing,
                                   lavaan_missing = ts_missing,
                                   test = attr(res, "model_test"),
                                   rows = ts_rows,
                                   footer = ts_footer)
  }
  attr(res, "model") <- syntax
  res
}


# The versions of the test and of the cfi and rmsea columns of a ladder,
# lavaan's names for each, most preferred first.
.MI_TEST_VERSIONS <- list(scaled = c("chisq.scaled", "df.scaled",
                                     "pvalue.scaled"),
                          standard = c("chisq", "df", "pvalue"))
.MI_INDEX_VERSIONS <- list(robust = c("cfi.robust", "rmsea.robust"),
                           scaled = c("cfi.scaled", "rmsea.scaled"),
                           standard = c("cfi", "rmsea"))

# The candidate versions of the test (`what` = "test") or of the indices
# (`what` = "indices") of a ladder, most preferred first, given whether its
# fit carries a scaled statistic and how its indices fall back (`method`,
# "fiml", "two_stage", "ml", or "other"; see
# .measurement_invariance_indices() for the order).
.measurement_invariance_candidates <- function(what, scaled, method) {
  if (what == "test") {
    return(if (scaled) .MI_TEST_VERSIONS else .MI_TEST_VERSIONS["standard"])
  }
  keep <- switch(method,
                 fiml = c("robust", if (scaled) "scaled", "standard"),
                 two_stage = c("robust", "scaled", "standard"),
                 ml = c(if (scaled) "robust", "standard"),
                 c(if (scaled) c("robust", "scaled"), "standard"))
  .MI_INDEX_VERSIONS[keep]
}

# The measures lavaan::fitMeasures() is asked for at every rung: every
# candidate of both kinds, and the baseline statistic, against which
# .sem_exact_fit() measures a statistic of zero. A name lavaan does not
# compute for a fit is left out of what it returns.
.measurement_invariance_measures <- function(scaled, method) {
  unique(c(unlist(
    c(.measurement_invariance_candidates("test", scaled, method),
      .measurement_invariance_candidates("indices", scaled, method)),
    use.names = FALSE), "baseline.chisq"))
}

# The first of `candidates` (a named list of lavaan measure names) that
# every rung reports. A rung with degrees of freedom must report it as a
# finite number. A saturated rung (df = 0 and a statistic of zero,
# .sem_exact_fit()) fits exactly in every version and does not decide:
# lavaan forms some of its indices there from a ratio of zeros (the robust
# CFI of a saturated model is 1 on one data set and NA on another, by
# rounding error) or leaves them out, and neither is a reason to change
# the version of every other rung, so the rung reports lavaan's value of
# the version chosen, or NA, and the caller gives its cfi and rmsea the
# values every version implies (1 and 0). A rung with no degrees of
# freedom whose statistic is not zero (.sem_zero_df_misfit(), as when an
# inequality bound binds) does not fit exactly. It decides the version of
# the indices through its CFI alone (`what` "indices"), since its RMSEA
# divides by the degrees of freedom and is undefined in every version, and
# it does not decide the version of the test (`what` "test"), since no
# version has a reference distribution on zero degrees of freedom (lavaan
# leaves the p-value of every version NA there). Through the first version
# of the rule for such a rung it decided both through every measure, so a
# robust ladder whose configural rung was one reported the standard
# statistic at every rung, under a model_test attribute of "standard",
# beside the scaled difference tests lavaan had computed. `fms` is the list
# of the rungs' fitMeasures(). Returns the name of the version.
.measurement_invariance_choose <- function(fms, candidates,
                                           what = c("indices", "test")) {
  what <- match.arg(what)
  exact <- vapply(fms, .sem_exact_fit, logical(1L))
  misfit <- vapply(fms, .sem_zero_df_misfit, logical(1L))
  decides <- !exact & !(misfit & what == "test")
  for (v in names(candidates)) {
    measures <- candidates[[v]]
    ok <- vapply(seq_along(fms), function(i) {
      fm <- fms[[i]]
      need <- if (misfit[i]) measures[1L] else measures
      !decides[i] ||
        (all(need %in% names(fm)) && all(is.finite(fm[need])))
    }, logical(1L))
    if (all(ok)) return(v)
  }
  names(candidates)[length(candidates)]
}

# The version of the chi_square, df, and p_chi_square columns of a ladder:
# the scaled statistic for a fit that carries one, when lavaan reports it
# at every rung, and otherwise the standard one (lavaan reports no scaled
# statistic beside the bootstrap test = "bollen.stine", for example).
# Returns the three measure names and the version, the model_test
# attribute of a ladder not estimated in two stages.
.measurement_invariance_test_columns <- function(fms, scaled) {
  candidates <- .measurement_invariance_candidates("test", scaled, "other")
  version <- .measurement_invariance_choose(fms, candidates, "test")
  list(measures = candidates[[version]], version = version)
}

# The version of the cfi and rmsea columns of a ladder, one version at
# every rung (.measurement_invariance_choose()), so that the differences
# between rungs compare like with like. A fit with a scaled statistic
# reports lavaan's robust indices, as fa_confirmatory() does for the same model, and
# a fit with neither a scaled statistic nor FIML the standard ones. Where
# lavaan does not compute the robust indices at every rung, every rung
# falls back to one other version, which depends on `method`:
#
#   "ml"        a maximum likelihood fit of continuous items with the cases
#               it reads complete (complete data, listwise or pairwise
#               deletion). Its robust indices estimate the population
#               values of the ML indices with their small sample bias
#               corrected for nonnormality (Brosseau-Liard, Savalei, & Li,
#               2012; Brosseau-Liard & Savalei, 2014), and the standard
#               indices estimate the same population values with the bias
#               corrected as if the data were normal, so the ladder falls
#               back to them. The scaled indices, which put the scaled
#               statistic into the standard formulas, estimate population
#               values that change with the kurtosis of the data (Savalei,
#               2018, whose Eqs. 3 and 5 give the population values of the
#               ML indices, against those of the scaled RMSEA after them,
#               and who advises abandoning the scaled ones). lavaan computes
#               no robust indices beside the older mean and variance
#               adjusted statistic (test = "mean.var.adjusted"), whose
#               robust indices DMAR computes (.sem_robust_mva(); the caller
#               adds them to `fms` and passes dmar_robust = TRUE), or a
#               bootstrap test.
#   "other"     a fit of ordered items, or by least squares, whose standard
#               statistic is not the basis of an index, so the ladder falls
#               back to the scaled indices of a scaled fit. lavaan 0.7-3
#               leaves the robust indices of the ordered ladders the tests
#               fit under the delta parameterization NA at every rung that
#               constrains the thresholds.
#   "fiml"      full information maximum likelihood. lavaan's robust
#               versions are the corrected indices of Zhang and Savalei
#               (2023) whatever the test: the complete data fit function
#               evaluated at the FIML estimates, with the saturated FIML
#               estimates in place of the sample moments, which approximate
#               the indices complete data would give (their version 6,
#               which assumes normality, or version 3, which does not, when
#               the fit has a scaled statistic, as with estimator = "MLR").
#               Neither the standard nor the scaled indices estimate the
#               complete data values, and the ladder falls back, as fa_confirmatory()
#               does for its rows, to the scaled versions of a scaled fit
#               or to the standard ones.
#   "two_stage" two-stage estimation, whose robust versions are the
#               corrected indices of the same authors (.ts_indices_label()),
#               falling back as under FIML.
#
# Returns the two measure names, the version ("robust", "scaled", or
# "standard"; the fit_indices attribute), and the line naming it: under
# FIML always, when the ladder falls back, and when its robust indices are
# DMAR's (NULL otherwise), which under two-stage estimation replaces the
# line .ts_footer() writes for the indices.
.measurement_invariance_indices <- function(fms, scaled,
                                            method = c("fiml", "two_stage",
                                                       "ml", "other"),
                                            dmar_robust = FALSE,
                                            allow_robust = TRUE) {
  method <- match.arg(method)
  candidates <- .measurement_invariance_candidates("indices", scaled, method)
  if (!allow_robust) candidates <- candidates[names(candidates) != "robust"]
  version <- .measurement_invariance_choose(fms, candidates)
  fallback <- version != names(candidates)[1L]
  which <- if (version == "scaled") {
    "put the scaled statistic into the standard formulas"
  } else switch(method,
    fiml = "change with the amount of missing data",
    two_stage = "treat the EM estimates as complete data",
    ml = paste0("estimate the same population values with the bias ",
                "corrected as if the data were normal (Savalei, 2018)"),
    "are computed from the standard statistic")
  line <- if (method == "fiml" && !fallback) {
    paste0("cfi, rmsea: the FIML corrected indices of Zhang and Savalei ",
           "(2023), version ",
           if (scaled) "3 (no normality assumption)" else "6 (normal theory)",
           ", which approximate the values complete data would give and ",
           "are not computed from the chi_square column; delta_cfi and ",
           "delta_rmsea compare them")
  } else if (fallback && method %in% c("fiml", "two_stage")) {
    paste0("cfi, rmsea: lavaan could not compute the ",
           if (method == "fiml") "FIML" else "two-stage",
           " corrected indices of Zhang and Savalei (2023) at every rung, ",
           "so these columns hold its ", version, " versions, which ", which)
  } else if (fallback) {
    paste0("cfi, rmsea: lavaan did not compute its robust versions at ",
           "every rung, so these columns hold its ", version, " versions, ",
           "which ", which, "; delta_cfi and delta_rmsea compare them")
  } else if (dmar_robust && version == "robust") {
    paste0(.sem_robust_mva_line("cfi, rmsea"),
           "; delta_cfi and delta_rmsea compare them")
  }
  list(measures = candidates[[version]], version = version,
       lines = if (!is.null(line)) c(fit_columns = line))
}


# The line under the table of a ULS ladder of continuous items whose
# indices are formed from the ULS statistic (see the fit index section of
# measurement_invariance()): `version` is the version of the cfi and rmsea
# columns, `restricted` the levels of the rungs that restrict the means,
# and `rmsea_withheld` whether any RMSEA was set to NA. NULL when nothing
# was. `version` "browne" is the standard version of a ladder tested by
# Browne's residual based statistic, whose `restricted` rungs are those
# that restrict the means through parameters the covariances share.
.measurement_invariance_uls_line <- function(version, restricted,
                                             rmsea_withheld) {
  if (!rmsea_withheld && !length(restricted)) return(NULL)
  if (identical(version, "browne")) {
    several <- length(restricted) > 1L
    return(paste0(
      "cfi, rmsea: formed from Browne's residual based statistic (the ",
      "chi_square column), which has no units, so they do not change when ",
      "every item is multiplied by the same number while the ULS estimates ",
      "do not, as at a rung that leaves the means free; at the ",
      paste(restricted, collapse = " and "),
      if (several) " rungs, which restrict" else " rung, which restricts",
      " the means through parameters the covariances share (the loadings, ",
      "through the latent means), the ULS fit function trades the squared ",
      "mean residuals against the squared covariance residuals, which are ",
      "in other units, so the estimates change with the units, and with ",
      "them the statistic and every index formed from it, and the cfi and ",
      "rmsea columns and their differences are NA there (the chi_square ",
      "column is kept, since its chi square reference distribution holds at ",
      "any consistent estimates (Browne, 1984, Proposition 4), though its ",
      "value changes with the units)"))
  }
  means <- function(what) {
    if (!length(restricted)) return("")
    paste0("; at the ", paste(restricted, collapse = " and "),
           if (length(restricted) > 1L) " rungs, which restrict" else
             " rung, which restricts",
           " the means, the ULS fit function adds the squared mean residuals ",
           "to the squared covariance residuals, which are in other units, ",
           "so the statistic, the loadings, and every index change with the ",
           "units, and ", what, " are NA there (the chi_square column and ",
           "the difference test are kept, since their reference distribution ",
           "does not depend on the units, though their values do)")
  }
  free <- if (length(restricted)) " while the means are free" else ""
  switch(version,
    robust = paste0("cfi, rmsea: ML-derived robust indices are withheld for ",
      "continuous ULS; a common-units check does not validate their ",
      "derivation (Savalei, 2018)"),
    scaled = paste0(
      "cfi, rmsea: the scaled indices, which do not change when every item ",
      "is multiplied by the same number", free,
      means("the cfi and rmsea columns and their differences")),
    standard = paste0(
      "cfi, rmsea: the standard statistic of a ULS ladder of continuous ",
      "items is (N - G) times the ULS fit function, in the units of the squared ",
      "covariances, and the CFI and RMSEA formed from it set its degrees of ",
      "freedom against it, so they change with the units of the items and ",
      "the columns are NA"))
}


# The line under the table of a ladder with rungs that have no degrees of
# freedom and a statistic that is not zero (`levels`, the names of those
# rungs; .sem_zero_df_misfit()), whose rmsea is NA.
.measurement_invariance_zero_df_line <- function(levels) {
  several <- length(levels) > 1L
  paste0("rmsea: NA at the ", paste(levels, collapse = " and "),
         if (several) " rungs, which have" else " rung, which has",
         " zero degrees of freedom and a statistic that is not zero (an ",
         "inequality constraint that binds costs no degree of freedom in ",
         "lavaan's count), so ", if (several) "they do" else "it does",
         " not fit exactly and the RMSEA, which divides by the degrees of ",
         "freedom, is undefined; the cfi is formed from the statistic, and ",
         "delta_rmsea at the next rung is NA")
}


# The lines under the table: those every table carries (the two-stage
# lines among them), and the line naming the version of the cfi and rmsea
# columns of a ladder fit by full information maximum likelihood or of a
# scaled ladder that fell back from the robust versions (the fit_columns
# attribute; .measurement_invariance_indices()).
#' @exportS3Method .dmar_tbl_footer dmar_measurement_invariance
#' @noRd
.dmar_tbl_footer.dmar_measurement_invariance <- function(x, ...) {
  .dmar_footer_add(NextMethod(), attr(x, "fit_columns", exact = TRUE),
                   "fit_columns")
}


# Observed categories per declared ordered indicator. The ladder depends on
# the count: a two-category item has one threshold, which is not separately
# identified from the intercept of its underlying response, and no free
# residual variance, so a thresholds rung and a strict rung would impose no
# restriction on it. Levels a factor carries but the data never realize are
# not categories, hence the count of distinct non-missing values rather than
# nlevels().
.measurement_invariance_n_categories <- function(data, items) {
  vapply(items, function(item) {
    x <- data[[item]]
    length(unique(x[!is.na(x)]))
  }, integer(1L))
}


# The missing data treatment of a ladder. The default of `missing` is
# "fiml" since DMAR 1.1.0 (DMAR 1.0.x deleted incomplete cases listwise),
# as in fa_confirmatory(), and it resolves by the rule of the reliability
# family (.relia_resolve_missing()). When no value of a model variable is
# missing among the cases with an observed group (lavaan sets the others
# aside under either treatment), the ladder is fitted as it stands
# ("listwise", which then deletes nothing), so that a complete data ladder
# is the one DMAR 1.0.x fitted: under FIML the ML rows would be the same,
# and the robust cfi and rmsea of an MLR ladder would be the corrected
# version 3 indices of Zhang and Savalei (2023), which differ from the
# complete data ones. A case with no observed model variable does not make
# the data incomplete, since listwise deletion and FIML both set it aside
# (lavaan drops it under FIML with a warning that it is empty). Incomplete
# data are fitted by FIML when the call allows it. FIML needs continuous items,
# estimator = "ML" or "MLR" (lavaan refuses MLM, MLMV, ULSM, and ULSMV
# under missing = "ml" and turns GLS, WLS, and ULS into two-stage
# estimation), standard errors other than robust.sem, and a test it can
# compute (lavaan falls back to listwise deletion for se = "robust.sem"
# and for the Satorra-Bentler kind of test). When the call asks for one of
# those, incomplete data stop it with the remedy named
# (.measurement_invariance_fiml_conflicts()), rather than deleting cases
# unasked. An explicit request is returned unchanged; an explicit FIML
# request that lavaan would not honor is refused by the caller
# (.measurement_invariance_check_fiml_request()).
#
# Inputs:
#   missing   : the validated value of `missing`.
#   default   : TRUE when the caller left `missing` at its default.
#   data      : the data frame.
#   variables : the observed variables of the model.
#   group     : the name of the grouping column.
#   ordered   : TRUE when ordered indicators were declared.
#   estimator : the estimator asked for.
#   se, test  : lavaan's se and test passed through ..., or NULL.
# Returns:
#   "listwise", "fiml", or `missing` itself.
.measurement_invariance_resolve_missing <- function(missing, default, data,
                                                    variables, group,
                                                    ordered, estimator,
                                                    se = NULL,
                                                    test = NULL) {
  if (!default) return(missing)
  present <- intersect(variables, names(data))
  block <- data[!is.na(data[[group]]), present, drop = FALSE]
  block <- block[rowSums(!is.na(block)) > 0L, , drop = FALSE]
  incomplete <- length(present) > 0L && anyNA(block)
  if (!incomplete) return("listwise")
  conflicts <- .measurement_invariance_fiml_conflicts(
    estimator, se, test, ordered = ordered, check_estimator = TRUE)
  if (!length(conflicts$why)) return("fiml")
  .dmar_refuse_argument(c(
    "Full information maximum likelihood, the default of 'missing', is ",
    "not available with ", paste(conflicts$why, collapse = " and "),
    ", and the items have missing values, so no cases are deleted ",
    "unasked. Use ",
    paste(c("missing = \"listwise\" to fit the complete cases",
            conflicts$remedy), collapse = ", or "),
    "."), arg = "missing")
}

# .measurement_invariance_fiml_conflicts(estimator, se, test,
#                                        ordered = FALSE,
#                                        check_estimator = FALSE)
#
# The options of a call that lavaan does not fit by full information
# maximum likelihood, and what to use in their place under FIML. lavaan
# falls back to listwise deletion for se = "robust.sem" and for a test of
# the Satorra-Bentler kind (.LAVAAN_LISTWISE_TESTS); ordered items take
# the categorical estimator, which has no FIML; and of the estimators only
# ML and MLR are fitted by FIML (check_estimator). A remedy for a test
# depends on the estimator: with MLR the scaled statistic FIML accepts is
# the estimator's own, so the remedy is to leave 'test' out; with ML or
# another estimator it is MLR in place of the test; with MLM or MLMV,
# whose own remedy is MLR, the two are named once.
#
# Inputs:
#   estimator       : the estimator asked for.
#   se, test        : lavaan's se and test passed through ..., or NULL.
#   ordered         : TRUE when ordered indicators were declared.
#   check_estimator : TRUE to count an estimator other than ML and MLR.
# Returns:
#   A list with `why` (the options, as the caller wrote them) and `remedy`
#   (phrases that follow "Use missing = \"listwise\" ...," in a message,
#   joined by ", or "; the one for MLR is a verb phrase of its own, "leave
#   out 'test'", and comes last); both empty when there is no conflict.
.measurement_invariance_fiml_conflicts <- function(estimator, se, test,
                                                   ordered = FALSE,
                                                   check_estimator = FALSE) {
  why <- character()
  if (ordered) why <- c(why, "ordered items")
  if (check_estimator && !estimator %in% c("ML", "MLR")) {
    why <- c(why, sprintf("estimator = \"%s\"", estimator))
  }
  robust_sem <- identical(se, "robust.sem")
  if (robust_sem) why <- c(why, "se = \"robust.sem\"")
  bad_test <- intersect(tolower(as.character(test)), .LAVAAN_LISTWISE_TESTS)
  if (length(bad_test)) {
    why <- c(why, sprintf("test = \"%s\"", bad_test))
  }
  mlm <- check_estimator && estimator %in% c("MLM", "MLMV")
  remedy <- c(
    if (ordered) "missing = \"pairwise\" for ordered items",
    if (mlm) {
      paste0("estimator = \"MLR\"", if (length(bad_test)) " without 'test'",
             ", whose robust standard errors and scaled test FIML accepts")
    },
    if (robust_sem) {
      paste0("se = \"robust.huber.white\" in place of se = ",
             "\"robust.sem\", which FIML accepts")
    },
    if (length(bad_test) && estimator == "MLR") {
      paste0("leave out 'test', since estimator = \"MLR\" computes its ",
             "own scaled statistic under FIML")
    } else if (length(bad_test) && !mlm) {
      paste0("estimator = \"MLR\" in place of the test, whose scaled ",
             "statistic FIML accepts")
    })
  list(why = why, remedy = remedy)
}

# .measurement_invariance_check_fiml_request(missing, default, estimator,
#                                            se, test)
#
# An explicit missing = "fiml" (or "ml") with se = "robust.sem" or a test
# of the Satorra-Bentler kind passed through ...: lavaan fits those by
# listwise deletion in place of FIML, with a warning, so the ladder would
# not be the one asked for. Refused before any fit, whether or not a value
# is missing, as fa_confirmatory() refuses it; through version 1.0.1 the
# ladder was the complete cases' one. Any other route by which
# lavaan replaces FIML (another spelling of such a test, or GLS, WLS, and
# ULS, which it fits in two stages) is caught on the fits themselves
# (.measurement_invariance_check_fiml_fits()).
#
# Inputs:
#   missing   : the resolved value of `missing`.
#   default   : TRUE when the caller left `missing` at its default (the
#               default refuses these options itself, and only when a
#               value is missing; .measurement_invariance_resolve_missing()).
#   estimator : the estimator asked for.
#   se, test  : lavaan's se and test passed through ..., or NULL.
# Returns:
#   NULL, invisibly; otherwise an error.
.measurement_invariance_check_fiml_request <- function(missing, default,
                                                       estimator, se,
                                                       test) {
  if (default || !missing %in% c("ml", "fiml")) return(invisible(NULL))
  conflicts <- .measurement_invariance_fiml_conflicts(estimator, se, test)
  if (!length(conflicts$why)) return(invisible(NULL))
  .dmar_refuse_argument(c(
    "missing = \"", missing, "\" asks for full information maximum ",
    "likelihood, which lavaan does not compute with ",
    paste(conflicts$why, collapse = " and "), ": it fits the complete ",
    "cases instead, whether or not a value is missing, so the ladder ",
    "would not be the one asked for. Use ",
    paste(c(paste0("missing = \"listwise\" to fit the complete cases ",
                   "with ",
                   if (length(conflicts$why) > 1L) "them" else "it"),
            conflicts$remedy),
          collapse = ", or "),
    "."), arg = NULL)
}

# .measurement_invariance_check_fiml_fits(fits, missing)
#
# The fits of a ladder that asked for FIML, held to it: lavaan replaces
# FIML with listwise deletion for some standard errors and tests (with a
# warning) and with two-stage estimation for GLS, WLS, and ULS (without
# one), and a ladder fitted so is not the one asked for. The known
# routes are refused before the fit
# (.measurement_invariance_check_fiml_request()); this catches the rest,
# as fa_confirmatory() checks its fit.
#
# Inputs:
#   fits    : the lavaan fits of the rungs.
#   missing : the resolved value of `missing`.
# Returns:
#   NULL, invisibly; otherwise an error.
.measurement_invariance_check_fiml_fits <- function(fits, missing) {
  if (!missing %in% c("ml", "fiml")) return(invisible(NULL))
  used <- vapply(fits, function(f) {
    as.character(lavaan::lavInspect(f, "options")$missing)[1L]
  }, character(1L))
  changed <- unique(used[!used %in% c("ml", "ml.x")])
  if (!length(changed)) return(invisible(NULL))
  .dmar_refuse_argument(c(
    "missing = \"", missing, "\" asks for full information maximum ",
    "likelihood, but lavaan fitted the ladder with missing = \"",
    changed[1L], "\" instead, for this combination of estimator, ",
    "standard errors, and test, so no ladder is returned from that fit. ",
    "Use estimator = \"ML\" or \"MLR\" with lavaan's default standard ",
    "errors and test, which are computed under FIML, or name the ",
    "treatment wanted (missing = \"listwise\" or \"two_stage\")."), arg = NULL)
}

# .measurement_invariance_check_coverage(data, variables, aux, group)
#
# Under FIML lavaan estimates, within each group, the means, variances, and
# covariances of the model's observed variables and of any auxiliary
# variables (for the saturated model the fit indices compare against, and
# for the ladder). A variable that no case of a group observes, or a pair
# that no case of a group observes together, leaves its moments in that
# group without a single observation, and lavaan stops with an internal
# error ("subscript out of bounds", after a warning that some pairs have
# zero coverage). holzinger_swineford is such a data set: its
# t25_paper_form_board_r and t26_flags were given in Grant-White only.
# Through version 1.0.1 an explicit FIML request on such data met that
# error, and with FIML the default since 1.1.0 it would be the ordinary
# path. The variables are refused by name before any fit, with what
# to do instead, in an error of class dmar_group_coverage (a dmar_refusal,
# .dmar_refuse()) that carries the variables and pairs at issue.
#
# Listwise deletion is never the remedy: it keeps only the cases observed
# on every variable of the model, and a group that never observes one, or
# never observes two together, has no such case, so lavaan would stop on
# that group instead. The remedies are to leave the variable (or one of
# the pair) out of the model, or out of 'aux', for this comparison, and,
# when at least two groups would remain, to leave the groups at issue out
# of the data.
#
# Inputs:
#   data      : the data frame.
#   variables : the observed variables of the model.
#   aux       : the auxiliary variables passed through ..., or NULL.
#   group     : the name of the grouping column.
# Returns:
#   NULL, invisibly, when every group observes every variable and every
#   pair; otherwise an error.
.measurement_invariance_check_coverage <- function(data, variables, aux,
                                                   group) {
  model_vars <- intersect(variables, names(data))
  aux <- setdiff(intersect(aux, names(data)), model_vars)
  vars <- c(model_vars, aux)
  if (length(vars) < 1L) return(invisible(NULL))
  g <- data[[group]]
  keep <- !is.na(g)
  groups <- if (is.factor(g)) {
    levels(droplevels(g[keep]))
  } else {
    sort(unique(as.character(g[keep])))
  }
  unobserved <- data.frame(group = character(), variable = character(),
                           n_group = integer(), stringsAsFactors = FALSE)
  pairs <- data.frame(group = character(), variable_1 = character(),
                      variable_2 = character(), stringsAsFactors = FALSE)
  for (lev in groups) {
    rows <- keep & as.character(g) == lev
    obs <- !is.na(as.matrix(data[rows, vars, drop = FALSE]))
    colnames(obs) <- vars
    none <- vars[colSums(obs) == 0L]
    if (length(none)) {
      unobserved <- rbind(unobserved, data.frame(
        group = lev, variable = none, n_group = sum(rows),
        stringsAsFactors = FALSE))
    }
    seen <- vars[!vars %in% none]
    if (length(seen) > 1L) {
      co <- crossprod(obs[, seen, drop = FALSE])
      zero <- which(co == 0 & upper.tri(co), arr.ind = TRUE)
      if (nrow(zero)) {
        zero <- zero[order(zero[, 1L], zero[, 2L]), , drop = FALSE]
        pairs <- rbind(pairs, data.frame(
          group = lev, variable_1 = seen[zero[, 1L]],
          variable_2 = seen[zero[, 2L]], stringsAsFactors = FALSE))
      }
    }
  }
  if (!nrow(unobserved) && !nrow(pairs)) return(invisible(NULL))
  rownames(unobserved) <- NULL
  rownames(pairs) <- NULL
  facts <- c(
    vapply(unique(unobserved$group), function(lev) {
      u <- unobserved[unobserved$group == lev, , drop = FALSE]
      paste0("group ", lev, " never observes ",
             .measurement_invariance_and(u$variable), " (none of its ",
             u$n_group[1L], " cases has a value)")
    }, character(1L)),
    vapply(unique(pairs$group), function(lev) {
      p <- pairs[pairs$group == lev, , drop = FALSE]
      paste0("no case of group ", lev, " observes ",
             paste0(p$variable_1, " and ", p$variable_2, collapse = ", or "),
             " together")
    }, character(1L)))
  # What to leave out: a model variable a group never observes, one model
  # variable of each pair of model variables, and an auxiliary variable
  # wherever one is at issue (the auxiliary of a pair that has one).
  is_model <- function(v) v %in% model_vars
  model_pair <- is_model(pairs$variable_1) & is_model(pairs$variable_2)
  leave_model <- unique(unobserved$variable[is_model(unobserved$variable)])
  leave_aux <- unique(c(unobserved$variable[!is_model(unobserved$variable)],
                        pairs$variable_1[!is_model(pairs$variable_1)],
                        pairs$variable_2[!is_model(pairs$variable_2)]))
  pair_words <- if (!any(model_pair)) {
    NULL
  } else if (sum(model_pair) == 1L) {
    paste0(pairs$variable_1[model_pair], " or ", pairs$variable_2[model_pair])
  } else {
    "one variable of each such pair"
  }
  leave_words <- if (length(leave_model) && length(pair_words)) {
    paste0(.measurement_invariance_and(leave_model), ", and ", pair_words,
           ",")
  } else if (length(leave_model)) {
    .measurement_invariance_and(leave_model)
  } else {
    pair_words
  }
  model_groups <- unique(c(
    unobserved$group[is_model(unobserved$variable)],
    pairs$group[model_pair]))
  bad_groups <- unique(c(unobserved$group, pairs$group))
  remedy <- c(
    if (length(leave_words)) {
      paste0("leave ", leave_words, " out of the model for this ",
             "comparison, so that the ladder compares the groups on what ",
             "every group observes")
    },
    if (length(leave_aux)) {
      paste0("leave ", .measurement_invariance_and(leave_aux),
             " out of 'aux'")
    },
    if (length(groups) - length(bad_groups) >= 2L) {
      paste0("leave group", if (length(bad_groups) > 1L) "s", " ",
             .measurement_invariance_and(bad_groups), " out of 'data', ",
             "so that the ladder compares the other groups on every ",
             "variable")
    })
  remedy[1L] <- paste0(toupper(substr(remedy[1L], 1L, 1L)),
                       substring(remedy[1L], 2L))
  .dmar_refuse(c(
    "Full information maximum likelihood estimates the means, variances, ",
    "and covariances of the model's variables within each group, so ",
    "lavaan cannot fit the ladder when a group gives one of them no ",
    "observation, as here: ", paste(facts, collapse = "; "), ". ",
    paste(remedy, collapse = "; or "), ".",
    if (length(model_groups)) {
      paste0(" Listwise deletion (missing = \"listwise\") is no remedy: ",
             "it keeps only the cases observed on every variable of the ",
             "model, and group", if (length(model_groups) > 1L) "s", " ",
             .measurement_invariance_and(model_groups), " ",
             if (length(model_groups) > 1L) "have" else "has",
             " none, so lavaan would stop there instead.")
    }),
    "group_coverage", unobserved = unobserved, pairs = pairs)
}

# The stop for a rung lavaan's optimizer did not converge on: an error of
# class dmar_not_converged, the class fa_confirmatory() gives a fit that
# did not converge (R/fa_confirmatory_internals.R), so that a program
# catches it by name. Such a rung has no estimates to report, and its fit
# measures are not available. The message names the rung; the condition
# carries it as `level`, beside the lavaan code that reproduces the rung's
# fit (`syntax`, and `lavaan_call`, written by .no_max_call_text() with
# my_data in place of the data). lavaan reports
# the failure for the fit as a whole ("the optimizer warns that a
# solution has NOT been found!"), not for a group, so no group is named.
.measurement_invariance_not_converged <- function(level, syntax,
                                                  lavaan_call) {
  .dmar_refuse(
    paste0(
      "measurement_invariance() stopped: lavaan's optimizer did not ",
      "converge on the ", level, " rung, so that rung has no estimates ",
      "to report and no result is returned. DMAR fits each rung with ",
      "lavaan::cfa(), the estimation engine; see lavaan ",
      "for more details. Arguments passed through ... reach every rung, ",
      "among them start values, and std.lv = TRUE, which sets the scale ",
      "of each factor by fixing its variance rather than its first ",
      "loading. The error carries the lavaan code of the rung as its ",
      "syntax and lavaan_call elements (see ?measurement_invariance)."),
    "not_converged", level = level, syntax = syntax,
    lavaan_call = lavaan_call)
}

# "a", "a and b", "a, b, and c".
.measurement_invariance_and <- function(x) {
  x <- as.character(x)
  n <- length(x)
  if (n <= 1L) return(x)
  if (n == 2L) return(paste(x, collapse = " and "))
  paste0(paste(x[-n], collapse = ", "), ", and ", x[n])
}


# Builds lavaan measurement syntax from a named list of factors. Kept here
# rather than shared with fa_confirmatory() because fa_confirmatory() also has to encode its
# within-factor equality descriptors; this one only needs the pattern.
.measurement_invariance_syntax <- function(factors, data) {
  if (!length(factors) || is.null(names(factors)) ||
      any(names(factors) == "")) {
    .dmar_refuse_argument(c(
      "A list 'model' must be named: ",
      "list(factor_name = c(\"item_1\", \"item_2\", ...))."), arg = "model")
  }
  factor_names <- names(factors)
  if (anyDuplicated(factor_names)) {
    .dmar_refuse_argument("Factor names must be unique.", arg = NULL)
  }
  if (!all(factor_names == make.names(factor_names))) {
    .dmar_refuse_argument(c(
      "Factor names must be syntactically valid R names (they become ",
      "lavaan labels); rename, for example, 'my factor' to 'my_factor'."),
      arg = NULL)
  }
  for (f in factor_names) {
    if (!is.character(factors[[f]]) || length(factors[[f]]) < 2L) {
      .dmar_refuse_argument(
        sprintf("Factor '%s' must name two or more items.", f), arg = NULL)
    }
  }
  items <- unlist(factors, use.names = FALSE)
  if (anyDuplicated(items)) {
    .dmar_refuse_argument(c(
      "In the list form each item loads on exactly one factor (simple ",
      "structure); duplicated: ",
      paste(unique(items[duplicated(items)]), collapse = ", "),
      ". Give 'model' as lavaan syntax for cross-loadings."), arg = "model")
  }
  if (length(factors) == 1L && length(factors[[1L]]) < 3L) {
    .dmar_refuse_argument(
      "A single factor needs three or more items for identification.",
      arg = NULL)
  }
  # An item the data lack is refused as fa_confirmatory() refuses one its
  # factors name (.fa_items_absent(): "unknown", or "form" for an NA
  # name); the checks above concern the model rather than the naming of
  # its items and do not take that class.
  absent <- setdiff(items, names(data))
  if (length(absent)) {
    .fa_items_absent(c("Items not found in 'data': ",
                       paste(absent, collapse = ", "), "."), absent)
  }
  paste(vapply(factor_names, function(f)
    paste(f, "=~", paste(factors[[f]], collapse = " + ")), character(1L)),
    collapse = "\n")
}


# The corrected (scaled or mean and variance adjusted) test statistics a
# lavaan fit carries, given its test option. Since lavaan 0.7-3 an ML fit
# lists Browne's normal theory residual based statistic
# ("browne.residual.nt.model") beside the standard one by default, and a
# weighted MLR fit lists Browne's asymptotically distribution free residual
# statistic ("browne.residual.adf") ahead of the Yuan-Bentler one. Browne's
# statistics are alternative tests of the model computed from the
# residuals, not corrections of the standard statistic, so they are set
# aside with it; reading one as a scaled test would take an ML fit for a
# robust one. Returns character(0) for a fit with no corrected statistic.
.lavaan_corrected_tests <- function(test) {
  test <- as.character(test)
  test[!is.na(test) & !test %in% c("standard", "none") &
         !startsWith(test, "browne.residual")]
}


# Names the chi square difference test lavTestLRT() actually performed. The
# scaling method is chosen by lavaan to match the test in force, so it is
# read back from the returned object's heading rather than assumed.
.measurement_invariance_test_label <- function(lrt, scaled) {
  if (!scaled) return("standard chi square difference test")
  heading <- paste(attr(lrt, "heading"), collapse = " ")
  heading <- gsub("[\r\n]+", " ", heading)
  # lavaan's heading may quote the method with straight or curly quotes
  # depending on the locale and the fancy-quotes option, so both are matched.
  # The curly quotes are written as escapes because R code in a portable
  # package must be ASCII.
  pattern <- "method[[:space:]]*=[[:space:]]*[\"\u201c]([^\"\u201d]+)[\"\u201d]"
  method <- if (grepl(pattern, heading)) {
    sub(paste0(".*", pattern, ".*"), "\\1", heading)
  } else {
    NA_character_
  }
  if (is.na(method)) {
    "scaled chi square difference test"
  } else {
    paste0("scaled chi square difference test (method = \"", method, "\")")
  }
}
