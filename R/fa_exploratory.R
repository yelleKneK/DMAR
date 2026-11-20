#' Exploratory Common-Factor Analysis With Confidence Intervals
#'
#' Fits a common-factor model, with rotation and uncertainty for its rotated
#' loadings. This is an exploratory measurement model: every item can load
#' on every factor. The estimates and rotation-aware standard errors are
#' computed by \pkg{lavaan}, which is required; DMAR adds the table and
#' its intervals, the display, the refusals, and the fit measure policy.
#' The function does not compute principal components. It carries a name
#' of its own so that attaching \pkg{lavaan}, whose \code{efa()} it calls,
#' masks nothing. \code{efa_k()} is the same function under a shorter name,
#' as \code{cfa_k()} is for \code{\link{fa_confirmatory}}.
#'
#' @param data Numeric raw observations. Supply exactly one of \code{data}
#'   or \code{S}.
#' @param nfactors Number of common factors, a positive integer giving an
#'   identified exploratory model. With a \code{target}, its number of
#'   factors when \code{nfactors} is not supplied; a supplied value must
#'   agree with it.
#' @param S Named positive definite covariance matrix, with \code{N}. A
#'   correlation matrix is a covariance matrix with unit variances, and
#'   the original metric then coincides with the standardized one.
#' @param N Sample size for covariance input.
#' @param items The items to analyze, by name: columns of \code{data}, or
#'   rows and columns of \code{S}, which is then reduced to them, in the
#'   order given; only the reduced matrix need be finite and positive
#'   definite. \code{NULL}, the default, analyzes every column of
#'   \code{data} or of \code{S}, or, with a \code{target}, the items the
#'   target names.
#'   Names that cannot be read as at least three distinct items are refused
#'   with an error of class \code{dmar_items_invalid} (see
#'   \code{\link{dmar_metadata}}).
#' @param estimator \code{"ML"} (default) or \code{"MLR"} for continuous
#'   items; \code{"ULS"} for the minimum residual (minres) solution of
#'   Harman and Jones (1966), unweighted least squares fitted to the
#'   correlation matrix (see Details); explicitly request \code{"WLSMV"}
#'   for ordered items. MLR uses the same likelihood estimates with
#'   Huber-White sandwich errors.
#' @param missing \code{"fiml"} (default), \code{"ml"}, or
#'   \code{"listwise"} for continuous items under ML or MLR. With ordered
#'   items the default fits the items as they stand when no item value is
#'   missing, as \code{\link{fa_confirmatory}} does; with missing values
#'   ordered items require an explicit \code{"pairwise"} or
#'   \code{"listwise"} choice, and FIML is not substituted with deletion. ULS has no likelihood, so it
#'   takes complete raw data, an explicit \code{"listwise"} choice, or a
#'   correlation matrix estimated elsewhere through \code{S} (for example
#'   the maximum likelihood estimate of \code{\link{cov_fiml}}); with
#'   incomplete data the default \code{"fiml"} is refused rather than
#'   replaced. Covariance input has no missing-data route.
#' @param ordered \code{NULL}, \code{TRUE} for all items, or their names.
#' @param rotation The rotation criterion, by name. Oblique, the factors
#'   free to correlate: \code{"oblimin"} (default; direct oblimin with
#'   gamma 0, Jennrich & Sampson, 1966), \code{"quartimin"} (Carroll,
#'   1953), \code{"geomin"} (Yates, 1987; epsilon 0.001),
#'   \code{"simplimax"} (Kiers, 1994), \code{"oblimax"} (Saunders, 1961),
#'   \code{"infomax"} (McKeon, 1968, as Browne, 2001, presents it), and
#'   \code{"bentler"} (Bentler, 1977). Orthogonal, the factors kept
#'   uncorrelated: \code{"varimax"}
#'   (Kaiser, 1958; with Kaiser row weights, as lavaan applies it) and
#'   \code{"quartimax"} (Neuhaus & Wrigley, 1954). \code{"target"} and
#'   \code{"target_orthogonal"} rotate to the partially specified target
#'   given in \code{target}, the factors free to correlate (Browne, 1972b)
#'   or kept uncorrelated (Browne, 1972a); see Details. \code{"none"}
#'   returns the unrotated solution lavaan identifies, whose factors are
#'   uncorrelated. Each criterion is applied in the form named here by
#'   lavaan's gradient projection algorithm from 30 random starts; Browne
#'   (2001) reviews the family.
#' @param target For \code{rotation = "target"} or
#'   \code{"target_orthogonal"}, the partially specified target, in either
#'   of two forms (\code{NULL}, the default, for every other rotation). The
#'   usual one is a named list of item keys, the form
#'   \code{\link{fa_confirmatory}} takes for \code{factors}: each element
#'   names a factor and the items expected to load on it. An item keyed to
#'   a factor is targeted at 0 on every factor it is not keyed to and left
#'   unspecified on the factors it is keyed to (an item keyed to two
#'   factors is unspecified on both), and an analyzed item keyed to no
#'   factor is unspecified on every factor. The other is a numeric matrix
#'   of items by factors, with the items as row names in any order,
#'   \code{NA} for each unspecified loading, and elsewhere a target value
#'   for the standardized pattern loading, strictly between -1 and 1 (see
#'   Details); column names, if given, name the factors.
#' @param conf_level Confidence level of the intervals for the parameters,
#'   the communalities, and the determinacy and \emph{H}, default 0.95.
#' @param rmsea_conf_level Confidence level of the RMSEA interval on the
#'   fit line and in the \code{fit_measures} attribute (with the gamma hat
#'   and adjusted gamma hat intervals there, which come from the same
#'   noncentrality limits), default 0.90, the conventional level, whose
#'   limits decide the tests of close and not close fit at .05 (see
#'   \code{\link{sem_fit_measures}}); \code{rmsea_conf_level = 0.95} gives a
#'   95 percent interval. It is separate from \code{conf_level}, which does
#'   not reach the RMSEA.
#' @param seed Optional seed for rotation starts; default \code{NULL} uses
#'   the current random stream. A supplied seed restores the caller's stream.
#' @param point_estimates_only \code{FALSE} (default) for the full analysis;
#'   \code{TRUE} for the same estimates without their standard errors,
#'   tests, and intervals, which are then not computed (see Details), for
#'   a program that needs only the estimates and wants them sooner.
#' @param x For \code{print()}, a table returned by \code{fa_exploratory()}.
#' @param digits For \code{print()}, the decimal places of the display
#'   (default \code{getOption("dmar.digits", 3)}; see
#'   \code{\link{dmar_tbl}}).
#' @param \dots For \code{print()}, passed to \code{print.data.frame()}.
#'
#' @details
#' The model is fitted to the covariance matrix of the items, with their
#' means under FIML, and the standardized solution is derived from it
#' (Cudeck, 1989, on why a correlation matrix is not analyzed as if it
#' were a covariance matrix); the table carries both metrics. Two routes
#' are correlation structures by construction: the minimum residual
#' solution, which is defined on the correlation matrix and is not
#' invariant to the scale of the items, and the WLSMV analysis of ordered
#' items, which fits their polychoric correlations and thresholds.
#'
#' FIML fits the raw incomplete observations directly, under a normal
#' working likelihood and missingness at random given the included items.
#' It does not treat an estimated unrestricted covariance matrix as if it
#' were an ordinary complete-data sample covariance. Factor signs and
#' order are arbitrary. Wald intervals describe the selected rotated
#' solution, conditional on the requested factor count; they do not include
#' uncertainty from choosing that count or rotation.
#'
#' \strong{Display.} The printed result, and \code{knit_print()} and
#' \code{\link{as_kable}} of it, show the solution rather than the
#' parameter list: one row per item with its standardized pattern loading
#' on each factor and that loading's interval, its communality with its
#' interval, and its complexity; above the table the sample, estimator,
#' missing data method, rotation and its form, and standard error method;
#' beneath it the variance each factor accounts for, the factor
#' correlations with their intervals (or the statement that the rotation
#' fixed them at zero), each factor's determinacy and coefficient \emph{H}
#' with their intervals, and the fit. The table itself keeps every
#' parameter in both metrics, which \code{tidy()} returns, and a subset of
#' it prints as an ordinary DMAR table.
#'
#' Hofmann's (1978) complexity of an item is the square of the sum of its
#' squared pattern loadings divided by the sum of its loadings to the
#' fourth power: 1 for an item that loads on one factor alone, and the
#' number of factors for an item that loads equally on all of them. The
#' variance a factor accounts for is its diagonal element of
#' Phi Lambda' Lambda, with Lambda the standardized pattern loadings and
#' Phi the factor correlation matrix (the column sum of squared loadings
#' when the factors are uncorrelated); these sum to the sum of the
#' communalities, and each is also given as a proportion of the total
#' variance of the items (their number, in the standardized metric) and of
#' the common variance. No cumulative column is shown, since the order of
#' the factors is arbitrary. The \code{variance_accounted} attribute holds
#' this table. The fit line reads the \code{fit_measures} attribute: the
#' chi square statistic with its degrees of freedom and p-value, the RMSEA
#' with its interval at \code{rmsea_conf_level} (90 percent by default),
#' which the line names, the CFI, TLI, and SRMR, and the version lines that
#' say which indices they are (see \code{\link{sem_fit_measures}}); under
#' ULS it gives the discrepancy statistic, the NFI, and the SRMR (below).
#' The measures are computed from the unrotated solution of the same
#' model: the chi square statistic and the standard indices do not depend
#' on the rotation, but lavaan's FIML-corrected indices of a rotated object
#' do, since its rotated parameterization counts more free parameters than
#' the model has, so the fit of a rotated solution is taken from the
#' unrotated fit, and every rotation reports the same fit.
#'
#' \strong{Rotation to a partially specified target.} A partially
#' specified target states what is expected of some loadings and nothing
#' of the rest: each loading is either given a target value or left
#' unspecified. The rotation minimizes the sum of squared differences
#' between the specified loadings and their target values, in the
#' standardized metric, and asks nothing of the unspecified ones (Browne,
#' 2001, p. 124). Zeros are the usual specification, since which loadings
#' should be near zero is easier to say than how large the others should
#' be (Browne, 2001, p. 124); the list form of \code{target} specifies
#' exactly those, a zero for each item on each factor it is not keyed to.
#' Unlike a confirmatory model, which forces the specified loadings to
#' their values, the rotation brings them only as close to their targets
#' as the data allow, so a specified zero whose rotated loading is far
#' from zero signals a misspecified zero (Browne, 2001, pp. 124--125), and
#' its interval in the table shows by how much. A specified value must lie
#' strictly between -1 and 1, and a target with a value outside that range
#' is refused: the target is stated in the standardized metric, and such a
#' value is more likely a loading taken from an item's original metric
#' than an expected standardized loading. No row weights are applied
#' (Browne, 2001, p. 125, finds the Kaiser and Cureton-Mulaik row
#' weights unsuited to a partially specified target), and the rotation is
#' iterated to a gradient tolerance of 1e-8 rather than lavaan's 1e-5, so
#' that zeros the rotation attains exactly (below) are attained to about
#' that precision. The factors keep lavaan's names (f1, f2, ...) in the
#' order of the target's columns, and the display prints the target
#' beneath the table with the names it gave. Rotation does not change the
#' model, so the fit is that of every other rotation.
#'
#' When every specified value is zero, the sign of each factor is
#' arbitrary, and each factor is reflected, as lavaan does by default, so
#' that its loadings have a positive sum. Since lavaan reflects all the
#' factors or none, no factor is reflected when any specified value is
#' nonzero. A factor whose column holds a nonzero value then takes its
#' sign from the target (a negative target keeps its sign), but a factor
#' whose specified values are all zero takes the sign of the random start
#' from which the rotation reached its minimum. A supplied \code{seed}
#' fixes that sign; with \code{seed = NULL} the loadings of such a factor,
#' and its correlations with the other factors, can change sign from one
#' run to the next. Either sign gives the same solution, so such a factor
#' may be read with its loadings and its correlations reversed together.
#'
#' The target must identify the rotation, and lavaan does not check that
#' it does: an unidentified rotation comes back with finite standard
#' errors and no warning. DMAR checks the conditions of Algina (1980)
#' before reporting. In oblique rotation each factor's column needs at
#' least m - 1 specified zeros, for m factors, and the standardized
#' unrotated loadings of the items specified in it must have rank m - 1
#' (Browne, 2001, p. 115; Algina, 1980, Theorem 2). A column with a
#' nonzero target value needs one entry more, m entries with rank m
#' (Theorem 1): with the factor variances at one, m - 1 restrictions that
#' include a nonzero value need not identify the rotation, as Algina
#' (1980, pp. 393--394) shows by example, so DMAR asks for the larger
#' count. In orthogonal rotation fewer suffice: taken in increasing order
#' of their specified entries, the columns need at least 0, 1, ..., m - 1
#' zeros, and again a column with a nonzero value one entry more (Theorem
#' 3), and the rank of a column's rows is checked against the same
#' numbers. Two columns with the same specified entries and values would
#' leave their factors indistinguishable. A target that fails is refused
#' with an error of class \code{dmar_target_unidentified}, and one that
#' cannot be read as a target with class \code{dmar_target_invalid}; both
#' inherit \code{dmar_refusal} and carry a \code{reason} element (see
#' \code{\link{dmar_metadata}}). Only a set of rows that is rank deficient
#' to within the precision of the fitted loadings is refused (a singular
#' value below 1e-6 times the largest); one that is nearly deficient
#' shows itself in wide intervals.
#'
#' When a column carries only the minimum number of zeros, the rotation
#' can attain them exactly. Those loadings are then fixed by the target
#' rather than estimated (lavaan reports a standard error of zero and a
#' statistic that means nothing). The table shows each at its target value
#' itself, in place of the rotated estimate, which differs from it only by
#' the rotation's convergence error, with \code{se}, \code{statistic},
#' and \code{p_value} \code{NA} and both limits at the target, in both
#' metrics. When every specified loading is attained this way, as in an
#' oblique rotation with m - 1 zeros in every column, the rotated solution
#' is the confirmatory factor model with those loadings fixed at their
#' targets and the factor variances at one (the factors uncorrelated in an
#' orthogonal rotation).
#'
#' Myers, Jin, Ahn, Celimli, and Zopluoglu (2015) studied target rotation
#' by simulation. Three factors were measured by 10 or 20 items (p. 497).
#' The communalities were low (.20 to .40), high (.60 to .80), or wide
#' (.20 to .80), and minor factors added error to the model (both p. 498).
#' The target values were in error by up to two standard errors (p. 500).
#' The errors in the target values had a negligible effect on the accuracy
#' and the stability of the rotated loadings (pp. 502--503). Specifying
#' more targets per factor improved both somewhat, the accuracy most when
#' the communalities were low and the stability most when they were high
#' (p. 502).
#'
#' In their comparison with geomin, every target value was in error by two
#' standard errors (pp. 500--501). Target rotation recovered the
#' population loadings more accurately in all 27 comparisons, one for each
#' number of targets at each level of communality, though by small amounts
#' (p. 502). Their measure of accuracy was the absolute value of the bias
#' of a replication, the signed error averaged over all the loadings
#' (p. 501); its mean was at or below 0.10 for both rotations, and the
#' largest difference between them was 0.08 (p. 502). Geomin gave the
#' more stable solution in 20 of the 27, among them every comparison with
#' low communalities; target rotation was the more stable with high
#' communalities and six or more targets per factor, and with wide
#' communalities and nine or more (pp. 499, 502--503). Their advice on
#' the number of targets is to specify the target as fully as the theory
#' allows, particularly when the items are not uniformly highly reliable
#' indicators of the factors (p. 503). Comparing the rotations, they
#' conclude that the extra work of specifying a target may be worthwhile
#' for accuracy; for stability, they expect geomin to do better when the
#' items are less than uniformly reliable indicators, and target rotation
#' to do better when the items are uniformly reliable and several targets
#' are specified (p. 504). Less stability is more sampling variability,
#' which the intervals of a given analysis show.
#'
#' The intervals condition on a target fixed before the data were seen.
#' Browne (2001, pp. 125, 148) describes revising a target in sequence,
#' leaving unspecified the entries that proved misspecified and specifying
#' as zeros the loadings that proved near zero; such a sequence is
#' exploratory, capitalizes on chance, and calls for a confirmatory
#' analysis of new data. The \code{method} attribute holds the target
#' (\code{target}, items by factors, \code{NA} where unspecified, with the
#' target's factor names) and the number of loadings it specifies
#' (\code{n_specified}).
#'
#' \strong{Determinacy and H.} The determinacy of a factor is the
#' correlation of its regression factor score with the factor, the
#' multiple correlation of the factor with the items (Grice, 2001; the
#' factor determinacy index of Ferrando & Lorenzo-Seva, 2018):
#' \eqn{\sqrt{[\Phi \Lambda' R^{-1} \Lambda \Phi]_{kk}}}, with
#' \eqn{\Lambda} the standardized pattern loadings, \eqn{\Phi} the factor
#' correlations, and \eqn{R} the item correlation matrix the model
#' implies. Its square is the largest reliability any weighted sum of the
#' items reaches as a measure of the factor, which Ferrando and
#' Lorenzo-Seva call the generalized \emph{H} and
#' \code{\link{bifactor_indices}} the maximal reliability. Coefficient
#' \emph{H} is Hancock and Mueller's (2001) closed form on the factor's
#' standardized loadings \eqn{l_j} alone,
#' \eqn{H = (1 + (\sum_j l_j^2 / (1 - l_j^2))^{-1})^{-1}}, applied to
#' one factor of a multidimensional solution as Rodriguez, Reise, and
#' Haviland (2016) apply it, so that each item's variance from the other
#' factors counts as error. With one factor \emph{H} is the squared
#' determinacy; otherwise the two can differ in either direction, and the
#' squared determinacy is the quantity to read beside \emph{H}. Both are
#' reported for every factor with a delta method standard error over
#' lavaan's covariance matrix of the standardized estimates (the matrix
#' behind the standardized loading and communality intervals) and Wald
#' limits formed on the logit scale and transformed back (Browne, 1982;
#' the confidence interval transformation principle of Kelley, 2007: the
#' limits of a monotone transformation of a parameter are the
#' transformation of its limits), as \code{\link{reliability_H}} forms them
#' by default, so they lie within (0, 1); the \code{se} column is the standard error of the index
#' itself, and with one factor it is the standard error
#' \code{\link{fa_confirmatory}} gives its \emph{H} row for the same model
#' fitted the same way (by FIML). For ordered items both describe the
#' latent responses behind the categories, not composites of the observed
#' scores (the G-H-latent of Ferrando and Lorenzo-Seva), and the display
#' says so.
#' Beside each determinacy \eqn{\rho} the display and the
#' \code{factor_quality} attribute give \eqn{2\rho^2 - 1}
#' (\code{min_score_correlation}, with its delta method standard error and
#' the determinacy's limits mapped through it), the lowest correlation two
#' sets of factor scores equally consistent with the model can have
#' (Guttman, 1955, as Grice, 2001, p. 435, presents it): zero when
#' \eqn{\rho^2 = .5}, negative below. No line is drawn.
#' \code{\link{factor_scores}} gives the determinacy of the Bartlett score
#' as well, and for each missing data pattern. The \code{factor_quality}
#' attribute holds this table.
#'
#' \code{estimator = "ULS"} is the minimum residual method: the loadings
#' and unique variances minimize the sum of squared residuals of the
#' correlation matrix, the unique variances absorbing the diagonal, so that
#' the off-diagonal correlations are reproduced as closely as the common
#' factors allow (Harman & Jones, 1966; \code{psych::fa(fm = "minres")}
#' gives the same loadings and communalities, up to the order and signs of
#' the factors). The estimates need no
#' distributional assumption, and the correlation matrix is the structure
#' fitted (lavaan's \code{correlation = TRUE}), so rescaling an item does
#' not change them. With raw data the standard errors are lavaan's sandwich
#' for a least squares estimator with the normal theory covariance matrix
#' of the correlations (its \code{se = "robust.sem.nt"}; Browne, 1984),
#' and the rotation enters them as for ML. A supplied correlation or
#' covariance matrix carries no such covariance matrix, so with \code{S} the
#' ULS estimates and communalities are returned with their standard errors
#' and limits \code{NA}, and the note under the table says so. The ULS
#' statistic is in the units of the squared correlations and has no chi
#' square reference distribution, so in the \code{fit_measures} attribute
#' the \code{chi_square} row holds N - 1 times the ULS discrepancy with no
#' p-value, the RMSEA rows and the CFI, TLI, RNI, and IFI are \code{NA},
#' and the NFI and the SRMR, which do not rest on that distribution, are
#' reported (see \code{\link{sem_fit_measures}}); a maximum likelihood fit
#' of the same items gives the fit indices and their inference.
#'
#' The table contains original-unit parameters, standardized loadings and
#' factor correlations, and communalities with intervals. A communality is
#' one minus the standardized residual variance, so its interval reverses
#' the limits of that variance. Wald limits are not clipped to [0, 1].
#' Improper or unidentified solutions are refused, with an error of class
#' \code{dmar_efa_no_solution} whose \code{reason} element says whether
#' the model did not converge, converged to an improper solution (lavaan's
#' post-estimation check failed: a residual variance below zero, or a
#' covariance matrix that is not positive definite), or left the
#' rotation-aware standard errors unavailable; \code{\link{factor_count}}
#' reads it to report such a count as a row. When a residual variance is
#' below zero and the likelihood has no maximum at all (the likelihood
#' flat or still rising as the variance falls, as when too many factors
#' are extracted), the refusal says so instead (\code{reason}
#' \code{"no maximum"}, and the class \code{dmar_no_maximum} as well; see
#' \code{\link{fa_bifactor}} for the check) and prints the
#' \code{lavaan::efa()} code of the unrotated fit, which selects
#' \code{items} from the data or covariance matrix given, where
#' \code{bounds = "pos.var"}, \pkg{lavaan}'s own default, or fewer factors
#' are the choices. When the unrestricted
#' FIML likelihood diagnostic fails but the constrained model and its
#' H0-based errors remain eligible, parameter results are qualified and
#' H1-dependent fit indices are withheld. Absence of a diagnosed obstacle
#' does not certify a unique global maximum.
#'
#' \strong{Point estimates only.} With \code{point_estimates_only = TRUE}
#' the table has the rows and the estimates of the full analysis (the
#' loadings, factor correlations, and residual variances in both metrics,
#' the communalities, Hofmann's complexity, and each factor's determinacy
#' and \emph{H}), and so do the \code{variance_accounted},
#' \code{factor_quality}, and \code{fit_measures} attributes, but no
#' standard error, test of a parameter, or interval is computed: the
#' \code{se}, \code{statistic}, \code{p_value}, \code{ci_lower}, and
#' \code{ci_upper} columns are \code{NA}, as are the standard errors and
#' limits of \code{factor_quality}. lavaan is asked for no standard errors
#' (its \code{se = "none"}), so neither the information matrix behind them
#' (with its sandwich under MLR, WLSMV, and ULS) nor the rotation-aware
#' errors of the rotated loadings are formed, nor the delta method errors
#' of the standardized solution, the communalities, the determinacy, and
#' \emph{H}; and the rotated fit is made without a test statistic, which
#' DMAR never reads, since the fit measures come from the unrotated fit.
#' That fit keeps its test, so every fit measure is the one the full
#' analysis reports, the robust statistic and the FIML-corrected indices
#' under MLR and the scaled and shifted statistic under WLSMV included.
#' \code{conf_level} then sets nothing, since the only interval left is the
#' RMSEA interval of the fit line, at \code{rmsea_conf_level}, and the
#' display states no other level. A loading of a target rotation that is attained
#' exactly is recognized by lying within 1e-6 of its target, in place of
#' its standard error of zero, and shown at the target as in the full
#' analysis; in rare cases (a chance of the order of 1e-5 for each
#' specified loading) an estimated loading falls that close to its target
#' by chance and is shown at it, which the full analysis would estimate,
#' within 1e-6 of the same value. The input and target checks, and the
#' refusals of a solution that did not converge or is improper and of a
#' likelihood without a maximum (whose code reproduces the full fit), are
#' made as before. The refusal of a solution whose rotation-aware standard errors
#' are unavailable is the full analysis's guard against a solution that is
#' not identified at its estimates, whose estimates are then not unique; it
#' is kept by inverting the expected information of the unrotated fit,
#' which is cheap, and, only when that fails, the information its standard
#' errors would use (the observed information under FIML); a solution
#' whose information cannot be inverted is refused with the same class and
#' reason. A failure confined to the
#' rotation's part of the rotation-aware errors is not detected, so in that
#' rare case the estimates are returned for a solution the full analysis
#' refuses. A FIML fit whose standard errors would need the unavailable
#' unrestricted (H1) information is not refused, since none are computed;
#' its fit measures are withheld as in the full analysis, with a warning of
#' the same class whose note speaks of the estimates alone. The display
#' says that standard errors and intervals were not computed, the result
#' carries the setting as its \code{point_estimates_only} attribute, and
#' the lavaan object in \code{fit} has no standard errors (and, when
#' rotated, no test).
#' The time saved depends on the estimator and the data. In one set of
#' timings with three factors under oblimin, the nine classic Holzinger and
#' Swineford tests took about a third less time under ML and MLR (0.79 s
#' against 0.51 s under ML) and half or more less under ULS and WLSMV
#' (0.70 s against 0.27 s for the tests cut into five ordered categories),
#' and the rotated solution of 18 items with 1,000 cases and 410 patterns
#' of missing values took 29 s against 17 s under MLR with FIML, and 18 s
#' against 9 s under ML. What remains is the estimation itself and the fit
#' measures, chiefly the robust statistic and the FIML-corrected indices of
#' the unrotated fit, which are kept; an unrotated analysis therefore saves
#' little or nothing (16.5 s against 16.7 s on those 18 items), and a
#' program that also fits the unrotated model (for its own fit measures,
#' say) saves on the rotated fit alone: the two together took 45 s against
#' 34 s.
#'
#' @return A \code{dmar_tbl} data frame (classes \code{dmar_fa_exploratory}
#'   and \code{dmar_factor_analysis}) with \code{term}, the parameter
#'   labels \code{lhs}, \code{op}, and \code{rhs}, \code{estimate},
#'   \code{se}, \code{statistic}, \code{p_value}, \code{ci_lower},
#'   \code{ci_upper}, \code{conf_level}, and \code{metric}: the parameters
#'   in the original metric, the standardized solution, the communalities
#'   (op \code{"communality"}), Hofmann's complexity (op
#'   \code{"complexity"}, with no interval), and each factor's determinacy
#'   and coefficient \emph{H} (ops \code{"determinacy"} and \code{"H"}, the
#'   factor as both \code{lhs} and \code{rhs}). Attributes retain \code{fit}
#'   (the lavaan object), \code{method} (sample size, estimator, standard
#'   error and missing data methods, the rotation and whether it is
#'   orthogonal, the fitted structure, and for a target rotation the target
#'   and the number of loadings it specifies), \code{variance_accounted},
#'   \code{factor_quality}, eligible \code{fit_measures} (a
#'   \code{\link{sem_fit_measures}} table, whose RMSEA interval is at
#'   \code{rmsea_conf_level}), \code{rmsea_conf_level}, and
#'   \code{point_estimates_only} (\code{TRUE} or \code{FALSE}, as called).
#'   \code{tidy()} preserves the parameter rows; \code{glance()} reports
#'   the analysis method and sample size, and both carry
#'   \code{point_estimates_only} as a column;
#'   \code{print()} shows the display described in Details.
#' @references
#' Algina, J. (1980). A note on identification in the oblique and
#'   orthogonal factor analysis models. \emph{Psychometrika, 45}(3),
#'   393--396. \doi{10.1007/BF02293911}
#'
#' Bentler, P. M. (1977). Factor simplicity index and transformations.
#'   \emph{Psychometrika, 42}(2), 277--295. \doi{10.1007/BF02294054}
#'
#' Browne, M. W. (1972a). Orthogonal rotation to a partially specified
#'   target. \emph{British Journal of Mathematical and Statistical
#'   Psychology, 25}(1), 115--120.
#'   \doi{10.1111/j.2044-8317.1972.tb00482.x}
#'
#' Browne, M. W. (1972b). Oblique rotation to a partially specified target.
#'   \emph{British Journal of Mathematical and Statistical Psychology,
#'   25}(2), 207--212. \doi{10.1111/j.2044-8317.1972.tb00492.x}
#'
#' Browne, M. W. (1982). Covariance structures. In D. M. Hawkins (Ed.),
#'   \emph{Topics in applied multivariate analysis} (pp. 72--141).
#'   Cambridge, UK: Cambridge University Press.
#'
#' Browne, M. W. (1984). Asymptotically distribution-free methods for the
#'   analysis of covariance structures. \emph{British Journal of
#'   Mathematical and Statistical Psychology, 37}(1), 62--83.
#'   \doi{10.1111/j.2044-8317.1984.tb00789.x}
#'
#' Browne, M. W. (2001). An overview of analytic rotation in exploratory
#'   factor analysis. \emph{Multivariate Behavioral Research, 36}(1),
#'   111--150. \doi{10.1207/S15327906MBR3601_05}
#'
#' Carroll, J. B. (1953). An analytical solution for approximating simple
#'   structure in factor analysis. \emph{Psychometrika, 18}(1), 23--38.
#'   \doi{10.1007/BF02289025}
#'
#' Cudeck, R. (1989). Analysis of correlation matrices using covariance
#'   structure models. \emph{Psychological Bulletin, 105}(2), 317--327.
#'   \doi{10.1037/0033-2909.105.2.317}
#'
#' Ferrando, P. J., & Lorenzo-Seva, U. (2018). Assessing the quality and
#'   appropriateness of factor solutions and factor score estimates in
#'   exploratory item factor analysis. \emph{Educational and Psychological
#'   Measurement, 78}(5), 762--780. \doi{10.1177/0013164417719308}
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
#' Harman, H. H., & Jones, W. H. (1966). Factor analysis by minimizing
#'   residuals (minres). \emph{Psychometrika, 31}(3), 351--368.
#'   \doi{10.1007/BF02289468}
#'
#' Hofmann, R. J. (1978). Complexity and simplicity as objective indices
#'   descriptive of factor solutions. \emph{Multivariate Behavioral
#'   Research, 13}(2), 247--250. \doi{10.1207/s15327906mbr1302_9}
#'
#' Jennrich, R. I. (1973). Standard errors for obliquely rotated factor
#'   loadings. \emph{Psychometrika, 38}(4), 593--604.
#'   \doi{10.1007/BF02291497}
#'
#' Jennrich, R. I., & Sampson, P. F. (1966). Rotation for simple loadings.
#'   \emph{Psychometrika, 31}(3), 313--323. \doi{10.1007/BF02289465}
#'
#' Kaiser, H. F. (1958). The varimax criterion for analytic rotation in
#'   factor analysis. \emph{Psychometrika, 23}(3), 187--200.
#'   \doi{10.1007/BF02289233}
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#'   Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kiers, H. A. L. (1994). Simplimax: Oblique rotation to an optimal
#'   target with simple structure. \emph{Psychometrika, 59}(4), 567--579.
#'   \doi{10.1007/BF02294392}
#'
#' McKeon, J. J. (1968). \emph{Rotation for maximum association between
#'   factors and tests}. Unpublished manuscript, Biometric Laboratory,
#'   George Washington University.
#'
#' Myers, N. D., Jin, Y., Ahn, S., Celimli, S., & Zopluoglu, C. (2015).
#'   Rotation to a partially specified target matrix in exploratory factor
#'   analysis in practice. \emph{Behavior Research Methods, 47}(2),
#'   494--505. \doi{10.3758/s13428-014-0486-7}
#'
#' Neuhaus, J. O., & Wrigley, C. (1954). The quartimax method: An analytic
#'   approach to orthogonal simple structure. \emph{British Journal of
#'   Statistical Psychology, 7}(2), 81--91.
#'   \doi{10.1111/j.2044-8317.1954.tb00147.x}
#'
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating
#'   bifactor models: Calculating and interpreting statistical indices.
#'   \emph{Psychological Methods, 21}(2), 137--150.
#'   \doi{10.1037/met0000045}
#'
#' Rosseel, Y. (2012). lavaan: An R package for structural equation
#'   modeling. \emph{Journal of Statistical Software, 48}(2), 1--36.
#'   \doi{10.18637/jss.v048.i02}
#'
#' Saunders, D. R. (1961). The rationale for an "oblimax" method of
#'   transformation in factor analysis. \emph{Psychometrika, 26}(3),
#'   317--324. \doi{10.1007/BF02289800}
#'
#' Yates, A. (1987). \emph{Multivariate exploratory data analysis: A
#'   perspective on exploratory factor analysis}. State University of New
#'   York Press.
#' @author Ken Kelley \email{kkelley@@nd.edu}
#' @examples
#' items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
#'            "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
#' two <- fa_exploratory(holzinger_swineford, nfactors = 2, items = items,
#'                       seed = 113)
#' two                       # the solution: loadings with intervals,
#'                           # communalities, complexity, variance, fit
#' attr(two, "variance_accounted")
#' attr(two, "factor_quality")  # determinacy and H with intervals
#' generics::tidy(two)[1:6, ]   # every parameter, both metrics
#'
#' # An orthogonal rotation, its estimates alone: the standard errors,
#' # tests, and intervals are not computed, which saves time.
#' fa_exploratory(holzinger_swineford, nfactors = 2, items = items,
#'                rotation = "varimax", seed = 113, point_estimates_only = TRUE)
#'
#' # Rotation to a partially specified target, through efa_k(), the same
#' # function under a shorter name, here of the minimum residual (minres)
#' # solution: unweighted least squares on the correlation matrix, with
#' # sandwich standard errors from the raw data. Each item is targeted at
#' # zero on the factor it is not keyed to, and the items and the number of
#' # factors are read from the target.
#' efa_k(holzinger_swineford, rotation = "target", estimator = "ULS",
#'       seed = 113, target = list(spatial = items[1:3], verbal = items[4:6]))
#' @seealso \code{\link{parallel_analysis}} for the factor count,
#'   \code{\link{fa_confirmatory}} for the confirmatory model,
#'   \code{\link{factor_scores}}, which takes the lavaan object in the
#'   \code{fit} attribute, and \code{\link{reliability_H}} for coefficient
#'   \emph{H} of one scale with bootstrap intervals.
#' @export
fa_exploratory <- function(data = NULL, nfactors = 1L, S = NULL,
                           N = NULL, items = NULL,
                           estimator = "ML", missing = "fiml",
                           ordered = NULL, rotation = "oblimin",
                           target = NULL,
                           conf_level = 0.95, rmsea_conf_level = 0.90,
                           seed = NULL, point_estimates_only = FALSE) {
  # Read before the argument `missing` is matched below, whose value then
  # shares the name of the base function.
  missing_default <- base::missing(missing)
  nfactors_supplied <- !base::missing(nfactors)
  .fa_require_lavaan()
  .fa_level(conf_level)
  .fa_level(rmsea_conf_level, "rmsea_conf_level")
  if (!is.logical(point_estimates_only) || length(point_estimates_only) != 1L ||
      is.na(point_estimates_only)) {
    .dmar_refuse_argument(
      "'point_estimates_only' must be TRUE or FALSE.",
      arg = "point_estimates_only")
  }
  point_only <- point_estimates_only
  .dmar_local_seed(seed)
  if (is.null(data) == is.null(S)) {
    .dmar_refuse_argument(
      "Supply exactly one of 'data' or 'S'.", arg = c("data", "S"))
  }
  estimator <- match.arg(estimator, c("ML", "MLR", "ULS", "WLSMV"))
  rotations <- .efa_rotations()
  rotation <- match.arg(rotation, rotations$name)
  rot <- rotations[rotations$name == rotation, ]
  targeted <- identical(rot$lavaan, "pst")
  if (targeted && is.null(target)) {
    .efa_target_refuse(sprintf(paste0("rotation = \"%s\" rotates to the target given ",
                                      "in 'target', which is missing."), rotation),
                       "rotation")
  }
  if (!targeted && !is.null(target)) {
    .efa_target_refuse(paste0("A 'target' is used only by rotation = \"target\" ",
                              "(oblique) or \"target_orthogonal\"; choose one of ",
                              "them or drop the target."), "rotation")
  }
  # With no items named, the items a target names are analyzed, from raw
  # data or from S alike, as fa_confirmatory() analyzes the items its
  # factors name.
  target_items <- function(available, input) {
    named <- rownames(.efa_target_matrix(target))
    absent <- setdiff(named, available)
    if (length(absent)) {
      .efa_target_refuse(c("The target names items that are not ", input, ": ",
                           paste(absent, collapse = ", "), "."), "names")
    }
    named
  }
  if (targeted && !is.null(data) && is.null(items)) {
    items <- target_items(colnames(data), "columns of 'data'")
  }
  missing <- match.arg(missing, c("fiml", "ml", "listwise", "pairwise"))
  uls <- estimator == "ULS"
  ord <- character()
  # The variables of the input as given, before the items are selected: the
  # lavaan code of a refused fit selects the items from it when they are
  # not the whole of it (.efa_no_maximum()).
  input_names <- if (!is.null(data)) {
    colnames(as.data.frame(data))
  } else {
    colnames(S)
  }
  if (!is.null(data)) {
    data <- .fa_data(data, items)
    items <- names(data)
    if (isTRUE(ordered)) ord <- items else if (!is.null(ordered)) ord <- ordered
    if (!is.character(ord) || anyNA(ord) || any(!ord %in% items) ||
        anyDuplicated(ord)) .dmar_refuse_argument(
          "'ordered' must name analysis items.", arg = "ordered")
    # The missing data method, resolved and checked as factor_count()
    # resolves it (.efa_resolve_missing()): with ordered items and no value
    # missing, the default fits the items as they stand (listwise deletion
    # with nothing to delete); with a value missing, the choice between
    # pairwise and listwise deletion is the user's. ULS takes complete
    # data or an explicit listwise choice.
    missing <- .efa_resolve_missing(missing, missing_default, data, ord,
                                    estimator)
    if (missing == "listwise") data <- data[stats::complete.cases(data), , drop = FALSE]
    data <- data[rowSums(!is.na(data)) > 0L, , drop = FALSE]
    if (nrow(data) < 3L) .dmar_refuse_argument(
      "Fewer than three usable cases.", arg = NULL)
    args <- list(data = data, ordered = if (length(ord)) ord else NULL,
                 missing = if (uls) "listwise" else if (missing == "fiml") "ml" else missing)
    if (uls) args$correlation <- TRUE
  } else {
    # The items, named or (with a target) read from the names of S, reduce
    # S to their rows and columns before it is checked, as they select the
    # columns of raw data.
    if (targeted && is.null(items)) {
      items <- target_items(colnames(.fa_covariance_names(S)), "items of 'S'")
    }
    S <- .fa_covariance(S, N, items)
    if (!is.null(ordered) || !estimator %in% c("ML", "ULS")) {
      .dmar_refuse_argument(c(
        "Covariance-input EFA supports continuous ML or ULS only; robust ",
        "and ordered uncertainty requires raw data."), arg = NULL)
    }
    items <- colnames(S)
    if (uls) {
      # The correlation matrix is the structure ULS fits. A supplied matrix
      # carries no asymptotic covariance matrix for the sandwich standard
      # errors (lavaan needs the raw data or a NACOV), so they are not
      # computed and the table says so.
      args <- list(sample_cov = stats::cov2cor(S), sample_nobs = N,
                   correlation = TRUE, se = "none")
    } else {
      args <- list(sample_cov = S, sample_nobs = N)
    }
  }
  # With point_estimates_only, lavaan computes no standard errors: neither
  # the information matrix and its inverse (with the sandwich under MLR,
  # WLSMV, and ULS from raw data), nor the rotation-aware errors of a
  # rotated fit, nor the delta method errors of the standardized solution.
  # The estimates do not depend on them, and neither does any test
  # statistic or fit measure DMAR reports (lavaan computes the corrected
  # statistics of MLR and WLSMV from the information it needs whether or
  # not the errors are computed).
  if (point_only) args$se <- "none"
  p <- length(items)
  B <- NULL
  if (targeted) {
    # The target as a matrix of items by factors, NA where unspecified (see
    # R/fa_exploratory_internals.R). It gives the factor count when none is
    # supplied and must agree with one that is.
    B <- .efa_target_matrix(target, items)
    if (!nfactors_supplied) {
      nfactors <- ncol(B)
    } else if (!is.numeric(nfactors) || length(nfactors) != 1L ||
               !isTRUE(nfactors == ncol(B))) {
      .efa_target_refuse(sprintf(paste0("'nfactors' must equal the number of factors ",
                                        "in the target (%d), or be left out."), ncol(B)),
                         "nfactors")
    }
  }
  if (!is.numeric(nfactors) || length(nfactors) != 1L || !is.finite(nfactors) ||
      nfactors < 1L || nfactors >= p || nfactors %% 1 ||
      (p - nfactors)^2 < p + nfactors) {
    .dmar_refuse_argument(c(
      "'nfactors' must be a positive integer giving an identified EFA ",
      "model for these items."), arg = "nfactors")
  }
  efa_fit <- function(rotation_spec, extra = list()) {
    do.call(lavaan::efa, c(args, list(nfactors = as.integer(nfactors),
      estimator = estimator, rotation = rotation_spec, bounds = "none",
      output = "lavaan"), extra))
  }
  # An improper solution is refused (.efa_no_solution()), and first, when
  # an error variance is negative, the check whether the likelihood has a
  # maximum at all (.lavaan_no_maximum()): a fit on a ridge, where that
  # variance runs off below zero and the likelihood is flat or still
  # rising, is refused with the lavaan code that reproduces it
  # (.efa_no_maximum(), which selects the items from the input when they
  # are not the whole of it). The check reads the unrotated fit, whose
  # parameters are the ones lavaan estimated (a rotation changes neither
  # the likelihood nor the error variances), fitting it when the fit at
  # hand is rotated; it does not apply under ULS, which has no likelihood,
  # or WLSMV.
  select <- if (!identical(input_names, items)) items
  proper <- function(fit) {
    if (!isTRUE(lavaan::lavInspect(fit, "converged"))) {
      .efa_no_solution("not converged")
    }
    if (!isTRUE(lavaan::lavInspect(fit, "post.check"))) {
      unrotated <- if (identical(lavaan::lavInspect(fit, "options")$rotation,
                                 "none")) {
        fit
      } else {
        # The unrotated fit of the same model warns as the rotated one,
        # the fit the call asked for, already has: that the solution seems
        # to have run away, and that some estimated variances are
        # negative. Said twice they say nothing more, so they are muffled
        # at this one call.
        tryCatch(suppressWarnings(efa_fit("none")), error = function(e) NULL)
      }
      no_max <- if (!is.null(unrotated)) .lavaan_no_maximum(unrotated)
      if (!is.null(no_max)) {
        # The code reproduces the fit of the full analysis, with its
        # standard errors, whichever was asked for.
        .efa_no_maximum(no_max, c(args[names(args) != "se" | !point_only], list(
          nfactors = as.integer(nfactors), estimator = estimator,
          rotation = "none", bounds = "none", output = "lavaan")),
          select = select)
      }
      .efa_no_solution("improper")
    }
    fit
  }
  # The form of the rotation (orthogonal or oblique) is passed explicitly
  # and checked below, so a criterion is applied in the form the help
  # names and in no other.
  rotation_spec <- c(list(rot$lavaan),
                     if (rotation != "none") list(orthogonal = rot$orthogonal))
  fit_unrotated <- NULL
  if (targeted) {
    .efa_target_counts(B, rot$orthogonal)
    # The rank condition is read from the standardized unrotated loadings,
    # so the unrotated model, which the fit measures need in any case, is
    # fitted first (it draws no random numbers, so the rotation's starts
    # are those of any other order).
    fit_unrotated <- proper(efa_fit("none"))
    A <- lavaan::lavInspect(fit_unrotated, "std")$lambda
    if (!identical(rownames(A), items)) {
      .dmar_refuse_argument(
        "lavaan reordered the items. No substituted analysis is returned.",
        arg = NULL)
    }
    .efa_target_rank(B, A, rot$orthogonal)
    # The target and its mask are passed explicitly: lavaan's default mask
    # would leave every nonzero target value unspecified. Reflection, which
    # turns each column to a positive sum, is kept only when every
    # specified value is zero, where the sign of every factor is arbitrary.
    # lavaan reflects all the columns or none, so one nonzero value turns
    # it off for every factor: a column with a nonzero value takes its sign
    # from the target, and a column of zeros alone the sign of the random
    # start that reached the minimum (the help says so). No row weights
    # (Browne, 2001, p. 125, finds the Kaiser and Cureton-Mulaik row
    # weights unsuited to a partially specified target). The rotation is
    # iterated to a gradient tolerance of 1e-8 rather than lavaan's default
    # 1e-5, so that zeros a just identified target attains are attained to
    # about that precision (with the default, about 1e-5 in an orthogonal
    # rotation, and the solution then differs from the confirmatory model
    # it equals by about 1e-4).
    B0 <- B
    B0[is.na(B0)] <- 0
    mask <- 1 * !is.na(B)
    rotation_spec <- c(rotation_spec, list(target = unname(B0),
      target_mask = unname(mask), reflect = all(B0 == 0), row_weights = "none",
      gpa_tol = 1e-8))
  }
  # The test statistic of a rotated fit is never read: the fit measures
  # come from the unrotated fit (below). With point_estimates_only the
  # rotated fit is therefore made without a test as well. Under MLR with
  # FIML the robust statistic needs the observed information pattern by
  # pattern, and it was most of the rotated fit's time (at 18 items, 300
  # cases, and 165 missing data patterns, 4.7 s with errors and test, 2.9 s
  # without the errors, 0.3 s without either); the estimates are identical.
  rotated_extra <- if (point_only && rotation != "none") list(test = "none")
  fit <- proper(efa_fit(rotation_spec, rotated_extra))
  options <- lavaan::lavInspect(fit, "options")
  expected_se <- if (point_only) "none" else
    switch(estimator, ML = "standard", MLR = "robust.huber.white",
           ULS = if (is.null(data)) "none" else "robust.sem.nt",
           WLSMV = "robust.sem")
  if (!identical(options$se, expected_se) ||
      (!is.null(data) && !uls && missing %in% c("ml", "fiml") &&
         options$missing != "ml") ||
      (uls && (!isTRUE(options$correlation) ||
                 !identical(options$test, rotated_extra$test %||% "standard")))) {
    # The ULS story below (a statistic with no reference distribution, the
    # NFI and SRMR alone reported) holds for lavaan's standard statistic of
    # a correlation structure; a lavaan that tested it otherwise would need
    # the note and the tests revisited, not a silent change of report.
    .dmar_refuse_argument(c(
      "lavaan changed the requested EFA uncertainty, test, or missing-data ",
      "method. No substituted analysis is returned."), arg = NULL)
  }
  same_matrix <- function(a, b) {
    is.matrix(a) && isTRUE(all.equal(unname(a), unname(b), check.attributes = FALSE))
  }
  if (!identical(options$rotation, rot$lavaan) ||
      (rotation != "none" &&
         !identical(isTRUE(options$rotation.args$orthogonal), rot$orthogonal)) ||
      (targeted &&
         (!same_matrix(options$rotation.args$target, B0) ||
            !same_matrix(options$rotation.args$target_mask, mask) ||
            !identical(options$rotation.args$order_lv_by, "none") ||
            !identical(options$rotation.args$row_weights, "none")))) {
    # For a target rotation lavaan must also have kept the target, its
    # mask, and the order of the factors (its "pst" fixes the columns to
    # the target's, so f1 is the target's first factor), without row
    # weights.
    .dmar_refuse_argument(c(
      "lavaan changed the requested rotation or its form. No substituted ",
      "analysis is returned."), arg = NULL)
  }
  # The fit measures come from the unrotated solution of the same model:
  # the chi square statistic and the standard indices do not depend on the
  # rotation, but lavaan's FIML-corrected (robust) indices of a rotated
  # object do, since its rotated parameterization counts more free
  # parameters than the model has (on six Holzinger and Swineford items
  # the robust 95 percent RMSEA interval of one two-factor model ran
  # [0, 0.075] unrotated and [0, 0.154] under oblimin, with the same chi
  # square).
  if (is.null(fit_unrotated)) {
    fit_unrotated <- if (rotation == "none") fit else efa_fit("none")
  }
  if (uls && !identical(lavaan::lavInspect(fit_unrotated, "options")$test,
                        "standard")) {
    .dmar_refuse_argument(c(
      "lavaan changed the requested EFA uncertainty, test, or missing-data ",
      "method. No substituted analysis is returned."), arg = NULL)
  }
  diagnostic <- .fiml_fit_h1_diagnostic(fit)
  if (!is.null(diagnostic) && !.fiml_h0_vcov_supported(fit)) {
    .dmar_refuse(c("The requested EFA uncertainty requires unavailable H1 ",
                   "information. ", diagnostic),
                 "h1_unavailable", diagnostic = diagnostic)
  }
  # With point_estimates_only the check below of the rotation-aware
  # standard errors has nothing to read. Its role as the guard against a
  # solution that is not identified at its estimates is kept by inverting
  # the information of the unrotated fit, the cheap expected information
  # first, with the same class and reason (.efa_information_singular()).
  # ULS from a supplied matrix has no check in either mode.
  if (point_only && !(uls && is.null(data)) &&
      .efa_information_singular(fit_unrotated)) {
    .efa_no_solution("standard errors unavailable", point_only = TRUE)
  }
  pe <- lavaan::parameterEstimates(fit, level = conf_level)
  out <- .fa_parameter_table(pe, conf_level)
  out$metric <- "original"
  ss <- lavaan::standardizedSolution(fit, level = conf_level)
  ss$est <- ss$est.std
  standardized <- .fa_parameter_table(ss, conf_level)
  standardized$metric <- "standardized"
  no_se <- identical(expected_se, "none")
  if (!no_se && any(!is.finite(c(out$se[out$op == "=~"],
                                 standardized$se[standardized$op == "=~"])))) {
    .efa_no_solution("standard errors unavailable")
  }
  n_attained <- 0L
  if (targeted && (!no_se || (point_only && !(uls && is.null(data))))) {
    # A specified loading the rotation attains exactly is fixed by the
    # target, not estimated: lavaan gives it a standard error of zero, and
    # a statistic and p-value that mean nothing. It is shown at its target
    # value with both limits there and no standard error, statistic, or
    # p-value, in the standardized metric where the target is stated and,
    # for a zero, in the original metric as well (a nonzero standardized
    # value is a multiple of the item's model implied standard deviation
    # there, which is estimated). With point_estimates_only there is no
    # standard error to read, and a loading is taken as attained when it
    # lies within the rotation's convergence tolerance of its target
    # (.efa_target_attained()), so that the estimates are those of the
    # full analysis. ULS from a supplied matrix, which has no standard
    # errors either way, keeps the rotated estimates in both modes.
    factors <- lavaan::lavNames(fit, "lv")
    attained <- .efa_target_attained(standardized, B, factors, by_value = no_se)
    b <- attr(attained, "target")
    n_attained <- sum(attained)
    zero <- paste(standardized$lhs, standardized$rhs)[attained & b == 0]
    standardized <- .efa_fix_rows(standardized, attained, b[attained])
    out <- .efa_fix_rows(out, out$op == "=~" & paste(out$lhs, out$rhs) %in% zero, 0)
  }
  residual <- ss[ss$op == "~~" & ss$lhs == ss$rhs & ss$lhs %in% items, , drop = FALSE]
  communalities <- .fa_parameter_table(residual, conf_level)
  communalities$term <- paste(residual$lhs, "communality")
  communalities$op <- "communality"
  communalities$estimate <- 1 - residual$est
  residual_limits <- communalities[c("ci_lower", "ci_upper")]
  communalities$ci_lower <- 1 - residual_limits$ci_upper
  communalities$ci_upper <- 1 - residual_limits$ci_lower
  communalities$statistic <- communalities$p_value <- NA_real_
  communalities$metric <- "standardized"
  # Hofmann's complexity of each item and the variance each factor
  # accounts for, from the standardized pattern loadings and the factor
  # correlations (see R/fa_exploratory_internals.R).
  solution <- .efa_matrices(standardized)
  complexity <- communalities[match(solution$items, communalities$lhs), , drop = FALSE]
  complexity$term <- paste(solution$items, "complexity")
  complexity$lhs <- complexity$rhs <- solution$items
  complexity$op <- "complexity"
  complexity$estimate <- .efa_complexity(solution$L)
  complexity$se <- complexity$ci_lower <- complexity$ci_upper <- NA_real_
  # Each factor's determinacy and coefficient H, with delta method
  # standard errors and logit-scale limits (see .efa_factor_quality()).
  quality <- .efa_factor_quality(fit, ss, conf_level)
  out <- rbind(out, standardized, communalities, complexity, quality$rows)
  rownames(out) <- NULL
  if (no_se) {
    out$se <- out$statistic <- out$p_value <- NA_real_
    out$ci_lower <- out$ci_upper <- NA_real_
  }
  attr(out, "fit") <- fit
  attr(out, "method") <- c(list(N = sum(lavaan::lavInspect(fit, "nobs")),
    estimator = estimator, se = options$se,
    missing = if (is.null(data)) "covariance input" else
      if (options$missing %in% c("ml", "ml.x")) "fiml" else options$missing,
    rotation = rotation, orthogonal = rot$orthogonal,
    row_weights = options$rotation.args$row_weights, ordered = ord,
    structure = if (uls) "correlation" else "covariance"),
    if (targeted) list(target = B, n_specified = sum(!is.na(B))))
  attr(out, "variance_accounted") <- .efa_variance_accounted(solution$L, solution$Phi)
  attr(out, "factor_quality") <- quality$table
  if (is.null(diagnostic)) {
    # The ULS statistic has no reference distribution, which
    # sem_fit_measures() warns about and records under its table; here the
    # note below says so once, and the rows it withholds stay NA.
    attr(out, "fit_measures") <- withCallingHandlers(
      sem_fit_measures(fit_unrotated, conf_level = rmsea_conf_level),
      dmar_uls_standard_statistic = function(w) {
        if (uls) invokeRestart("muffleWarning")
      })
  }
  # The setting travels with the result, FALSE included, so that tidy(),
  # glance(), and a program reading the table can tell a table without
  # standard errors by request from one without them by necessity (ULS
  # from a supplied matrix).
  attr(out, "point_estimates_only") <- point_only
  # The level of the RMSEA interval, which differs from the table's
  # conf_level, so that dmar_metadata() and a program can read it beside
  # the level of the parameter intervals.
  attr(out, "rmsea_conf_level") <- rmsea_conf_level
  attained_text <- if (targeted && n_attained == sum(!is.na(B))) {
    paste0(" Every specified loading is attained exactly, so each is ",
           "fixed by the target, shown at its target value without a ",
           "standard error; the solution is the confirmatory model with ",
           "those loadings fixed.")
  } else if (targeted && n_attained) {
    sprintf(paste0(" %d of the %d specified loadings %s attained exactly ",
                   "and so fixed by the target, shown at the target value ",
                   "without a standard error."),
            n_attained, sum(!is.na(B)), if (n_attained == 1L) "is" else "are")
  }
  # The method is stated above the display (see
  # .dmar_tbl_footer.dmar_fa_exploratory()); the note under it carries
  # the qualifications alone.
  attr(out, "analysis_note") <- paste0(
    if (uls) {
      paste0("The minimum residual solution is unweighted least squares fitted ",
             "to the correlation matrix. The ULS statistic has no chi square ",
             "reference distribution: in the fit_measures attribute the ",
             "chi_square row is N - 1 times the ULS discrepancy without a ",
             "p-value, the RMSEA, CFI, TLI, RNI, and IFI are NA, and the NFI ",
             "and SRMR are reported. ")
    },
    if (point_only) {
      paste0("Standard errors, tests, and intervals were not computed ",
             "(point_estimates_only = TRUE); the estimates, the determinacy ",
             "and H, and the fit measures are those of the full analysis.",
             if (targeted && n_attained) {
               sprintf(paste0(" %d of the %d specified loadings %s attained ",
                              "exactly and shown at %s target value%s."),
                       n_attained, sum(!is.na(B)),
                       if (n_attained == 1L) "is" else "are",
                       if (n_attained == 1L) "its" else "their",
                       if (n_attained == 1L) "" else "s")
             })
    } else if (no_se) {
      paste0("Standard errors: none (a supplied matrix carries no covariance ",
             "matrix for the ULS sandwich; supply raw data for standard errors).")
    } else if (targeted) {
      paste0("Intervals condition on the selected factor count and on a target ",
             "fixed before the data were seen; a target revised after seeing a ",
             "rotation makes the analysis exploratory.", attained_text)
    } else {
      "Intervals condition on the selected factor count and rotation."
    })
  out <- .as_dmar_tbl(out, conf_level = conf_level, p_terms = "p_value",
                       subclass = "dmar_factor_analysis")
  class(out) <- c("dmar_fa_exploratory", class(out))
  .relia_h1_qualify(out, diagnostic, uncertainty = !point_only)
}

# The shorter name is the same function object, so the two cannot drift
# apart in their arguments or defaults.
#' @rdname fa_exploratory
#' @export
efa_k <- fa_exploratory

#' @rdname fa_exploratory
#' @export
print.dmar_fa_exploratory <- function(x, digits = getOption("dmar.digits", 3L), ...) {
  # A subset of the table (x[rows, cols]) keeps the class but not the
  # attributes the display reads, so it prints as an ordinary DMAR table.
  if (!.efa_display_ok(x)) return(invisible(NextMethod()))
  .dmar_tbl_write(x, digits = digits, print_args = c(list(right = TRUE), list(...)),
                  view_args = list(digits = digits))
}
