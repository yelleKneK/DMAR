#' Categorical Omega for Ordered-Categorical Items, With a Confidence Interval
#'
#' @description
#' Estimates categorical omega (\eqn{\omega_C}; Green & Yang, 2009;
#' Kelley & Pornprasertmanit, 2016) for a homogeneous composite of
#' ordered-categorical items and returns a bootstrap confidence interval.
#' The \code{denominator} argument selects whether the variance of the
#' summed categories is built from the unrestricted polychoric
#' correlations (\code{"observed"}, the default) or from the polychoric
#' correlations the fitted single-factor model reproduces
#' (\code{"model_implied"}, the coefficient as Green and Yang, 2009,
#' defined it). Incomplete rows are set aside by listwise deletion (the
#' default) or used through pairwise deletion or multiple imputation
#' (\code{missing}); under multiple imputation the interval is formed by
#' bootstrap then impute (Schomaker & Heumann, 2018).
#'
#' @details
#' Categorical omega is designed for items measured on an ordered
#' categorical scale (e.g., Likert items). It uses a probit-link
#' single-factor model in which each observed item \eqn{X_j} is modeled
#' as a categorization of an underlying continuous response variable
#' \eqn{X_j^{*}} via thresholds \eqn{t_{j,c}} (Muthén, 1984; Millsap &
#' Yun-Tein, 2004), fit by diagonally weighted least squares with mean-
#' and variance-adjusted test statistic (WLSMV), the standard estimator
#' for ordered categorical items. With the delta parameterization
#' (\eqn{Var(X_j^{*}) = 1}), the population categorical omega is
#' \deqn{\omega_C = \frac{\sum_{j=1}^{J} \sum_{j'=1}^{J} \sigma_{jj'}\!\left(\lambda_{j}\lambda_{j'}\right)}{\sum_{j=1}^{J} \sum_{j'=1}^{J} \sigma_{jj'}\!\left(\rho_{X_{j}^{*} X_{j'}^{*}}\right)},}
#' where \eqn{\sigma_{jj'}(r)} is the covariance of
#' \eqn{(X_{j}, X_{j'})} computed from a bivariate normal CDF over pairs
#' of category thresholds and a latent-response correlation \eqn{r}
#' (Green & Yang, 2009, Eq. 13--14, with the full coefficient their
#' Eq. 21; Kelley & Pornprasertmanit, 2016, Eq. 17--18). The numerator
#' uses the correlations of the common part, \eqn{\lambda_j
#' \lambda_{j'}}, under either setting of \code{denominator}; the two
#' settings differ in the correlations \eqn{\rho} used for the total
#' variance. The bivariate normal probabilities, one for each pair of
#' thresholds of each pair of items, are computed by \pkg{DMAR} itself,
#' by the method of Drezner and Wesolowsky (1990) as refined by Genz
#' (2004), which is the algorithm \pkg{mvtnorm} applies in two
#' dimensions. The two agree to within \eqn{4 \times 10^{-16}}, except
#' within about \eqn{10^{-7}} of a correlation of one, where
#' \pkg{mvtnorm} loses accuracy and \pkg{DMAR}'s computation does not.
#' Computed together, the probabilities take a small fraction of the time
#' of the factor model's fit.
#' \describe{
#'   \item{\code{"observed"} (default)}{The unrestricted polychoric
#'   correlations, estimated from the data by a saturated bivariate model
#'   (Kelley & Pornprasertmanit, 2016, Eq. 18). Structure that the
#'   single-factor model omits, such as a pair of items with shared
#'   wording, stays in the total variance, which makes this the
#'   categorical counterpart of robust omega (see
#'   \code{\link{reliability_omega}}). It is the coefficient
#'   \code{MBESS::ci.reliability(type = "categorical")} computes, apart
#'   from the samples described below in which the data do not determine
#'   a polychoric correlation.}
#'   \item{\code{"model_implied"}}{The latent-response correlations
#'   reproduced by the fitted single-factor model. This is categorical
#'   omega as Green and Yang (2009) defined it: their appendix program
#'   takes \dQuote{the reproduced (model-based) correlation matrix} as the
#'   input to the denominator, so correlated errors enter the total
#'   variance only when they are specified in the model.}
#' }
#' When the single-factor model reproduces the polychoric correlations
#' the two settings coincide, and they diverge as the model misfits.
#' Neither protects the numerator: shared structure among a few items
#' raises their loadings under either setting. The result therefore
#' reports the residual polychoric correlations (the unrestricted minus
#' the reproduced) and prints the largest under the table, so that a
#' high coefficient resting on redundant items is visible beside it.
#' When the shared structure sits on the items that dominate the factor,
#' the factor absorbs it and the residuals appear among the other items,
#' so the largest residual signals local dependence more reliably than it
#' locates it.
#'
#' \strong{What each denominator describes.} Both settings share the
#' numerator, the variance of the summed categories attributable to the
#' single factor, and differ in what they take the total variance of the
#' sum to be. The model implied setting takes the variance the fitted
#' single-factor model reproduces, so the coefficient describes the sum
#' score as the model describes it: covariance the model omits is absent
#' from the numerator and from the denominator alike. The observed
#' setting takes the variance implied by the polychoric correlations
#' estimated without the single-factor restriction (still under the
#' latent-response model that defines polychoric correlations), so
#' covariance the model omits is counted as variance of the sum that the
#' factor does not explain. The two estimate the same quantity when the
#' model holds. Kelley and Pornprasertmanit (2016) used the observed
#' setting because the total variance of the sum score is a property of
#' the data and need not be taken from a model that may be wrong; Green
#' and Yang (2009) used the model implied setting, which is the
#' reliability of the sum score under the fitted model.
#'
#' \strong{How the categories are scored.} Categorical omega is the
#' reliability of a sum in which each item's categories are scored as
#' consecutive ranks, 0, 1, \ldots, \eqn{C - 1}, in the order of the
#' item's codes: Green and Yang (2009, Eq. 11) and Kelley and
#' Pornprasertmanit (2016, Eq. 16) define the observed item as
#' \eqn{X_{ij} = c} for \eqn{c = 0, 1, \ldots, C - 1}, and \pkg{MBESS}
#' and \pkg{semTools} compute the same coefficient. The
#' function therefore uses the order of the codes in \code{data}, not
#' their values, and the coefficient is that of the sum of the codes as
#' given whenever every item's codes are equally spaced with one step
#' shared by all the items (codes 0, 2, 4, 6 on every item give the same
#' result as 1, 2, 3, 4). When they are not, as with codes -2, -1, 1, 2,
#' with items stepping by 2 beside items stepping by 1, or with a
#' category between two used ones that no respondent chose, the sum of
#' the codes is a different composite, and a message says so and names
#' the items; the coefficient is still that of the ranks. In the
#' examples behind this paragraph the coded sum's reliability was lower
#' by 0.009 to 0.012 for codes -2, -1, 1, 2 on every item and by 0.034
#' for two three-category items coded 1, 3, 5 beside two coded 1, 2, 3,
#' and an unused interior category changed it by about 0.001.
#'
#' \strong{When the data do not determine a polychoric correlation.} The
#' observed setting needs every unrestricted polychoric correlation.
#' When two items are skewed in opposite directions, as in Green and
#' Yang's (2009) mixed-skew conditions, their bivariate table can have
#' empty cells in a pattern (every respondent below the top category of
#' one item in the bottom category of the other) that leaves the
#' polychoric correlation without a finite standard error. The
#' unrestricted value is then arbitrary, and the single-factor fit gives
#' the correlation no weight, so the loadings rest on the correlations
#' the data do determine. Those identify the loadings only when they
#' connect every item to the others, directly or through other items,
#' and include at least one cycle of an odd number of items, such as
#' three items whose three correlations are all determined (the graph
#' whose edges are the determined pairs must be connected and not
#' bipartite); an item joined to the rest by a single determined
#' correlation is then identified as well. When every correlation between \code{y1} or \code{y2}
#' and \code{y3} or \code{y4} of the example below is undetermined, the
#' determined ones fix only the products of the loadings of
#' \code{y1} and \code{y2} and of \code{y3} and \code{y4}: the relative
#' sign of the two pairs is free, and \pkg{lavaan}'s default start
#' values can return loadings of opposite signs and a categorical omega
#' near 0.05, where the population value is 0.716. In 60 samples of 300
#' from that design, 13 had determined correlations that did not
#' identify the loadings: 2 with that sign reversal (0.05 and 0.09), 8
#' whose value moved with the start values by up to 0.06, and 3 whose
#' default start values gave an improper solution. The function
#' therefore checks the graph. When the determined
#' correlations identify the loadings, the observed setting stops with an
#' error naming the pair and pointing to the model implied setting, and
#' the model implied setting proceeds, says which pair it set aside, and
#' reports that pair's residual correlation as \code{NA}. When they do
#' not, both settings stop with an error naming the pairs and the blocks
#' of items no determined correlation joins, since the numerator takes
#' the loadings under either setting. The unrestricted polychoric
#' correlations come from the single-factor fit itself, which agrees with
#' a second, saturated fit whenever every correlation is determined.
#' \pkg{MBESS} computes them from such a saturated fit, which returns 0
#' for an undetermined correlation and so removes the pair from the total
#' variance without a message. The example below uses Green and Yang's
#' design at a sample size where samples rarely have this problem; at
#' \eqn{N = 500} most of them do.
#'
#' A bootstrap meets the problem far more often than the sample does.
#' The few respondents in the sparse cells of an oppositely skewed pair
#' are what determine its polychoric correlation, and a resample omits
#' some of them as a matter of course, so a sample whose correlations
#' are all determined still yields many resamples in which one is not.
#' On the data of the example below, which have no such pair, about
#' half of the replications of a bootstrap of the observed setting
#' have one (80, 95, and 97 of 200 in percentile bootstraps with
#' \code{seed} 1, 2, and 3). Those replications are set aside, counted in the
#' \code{bootstrap} attribute (\code{undetermined}; see \emph{Value}),
#' named in the printed line and in a warning, and the interval then
#' conditions on the replications that redrew the respondents who
#' determine the correlation, which is not the bootstrap distribution of
#' the coefficient; with more than half set aside, no interval is
#' reported (see below), and on these data the observed setting is
#' always near that point. The model implied setting needs no unrestricted
#' polychoric correlation and sets aside only the replications whose
#' determined correlations leave the loadings unidentified (counted as
#' \code{undetermined} as well), so with items skewed in opposite
#' directions it keeps more of them, but how many it keeps depends on
#' the sample size. In a percentile bootstrap of the example's complete
#' data (\eqn{N = 1000}, \code{B = 200}, \code{seed = 1}) it set aside 1
#' replication of 200, and over 20 samples of that size from the design
#' (\code{B = 50}) none in 15 and at most 18\%. Under pairwise deletion
#' of the example's incomplete copy it set aside 47 of 200: 27 whose
#' loadings were unidentified, 3 improper, and 17 whose fit failed, so
#' the caution applies. At \eqn{N = 300}, in the 45 samples of 60 whose
#' point estimate was identified and proper, it set aside a median of
#' 42\% of the replications (\code{B = 50}): 18\% where the sample
#' determined every correlation, 40\% where it left one undetermined,
#' and 64\% where it left two, and more than half in 15 of the 45.
#' When fewer than half of the replications return a value, no interval
#' is reported (see below). Under either denominator, a replication
#' whose determined correlations do not identify the loadings is set
#' aside rather than evaluated at arbitrary starting values.
#'
#' \strong{When the solution is improper.} With few respondents in the
#' extreme categories, or an item that nearly duplicates the factor,
#' the single-factor fit can return a standardized loading above one, a
#' negative residual variance (a Heywood case). Categorical omega is not
#' defined for such a solution, since its assembly evaluates the
#' bivariate normal distribution at the item's communality, and the
#' function stops with an error naming the item and its loading under
#' either setting of \code{denominator}. A bootstrap replication with an
#' improper solution is set aside and counted (\code{improper} in the
#' \code{bootstrap} attribute); these are the replications in which a
#' standardized loading ran past one (an error variance fell below
#' zero), so setting them aside can trim the upper tail of the bootstrap
#' distribution, which is why their number is reported.
#'
#' Kelley and Pornprasertmanit (2016) found in extensive Monte Carlo
#' simulation that the bias-corrected and accelerated (BCa) bootstrap
#' confidence interval for categorical omega achieved acceptable coverage
#' across a wide variety of threshold patterns, sample sizes, item
#' counts, and population reliability values, with one documented
#' exception: coverage dipped somewhat below the acceptable range when
#' the number of items and the population reliability were both high.
#' They specifically recommend BCa for categorical omega. Because no bootstrap runs in
#' \pkg{DMAR} unless the user requests one, the default output is the
#' point estimate with a message naming the call that produces the
#' recommended interval; request \code{ci_method = "bca"} to obtain it.
#' The interval is not a quick one. Every replication refits the
#' ordered-categorical factor model, whose sample statistics also
#' supply the unrestricted polychoric correlations behind the
#' denominator. That is a fraction of a second for a scale of moderate
#' length, so the recommended interval at the default
#' \code{B = 10000} can run for several minutes; when \code{B} is below
#' the number of cases, the jackknife behind the BCa acceleration adds
#' one such refit per case (at the default \code{B} the acceleration is
#' read from the replicates themselves at no extra cost; see the
#' \code{"bca"} entry below). The example on this page
#' therefore stops at the point estimate; the reported analysis adds
#' \code{ci_method = "bca"} to the same call, with \code{seed} for
#' reproducibility, and for these oppositely skewed items it bootstraps
#' the model implied setting, which at this sample size keeps nearly all
#' of the replications of the complete data (see above), or reports how
#' many replications of the observed setting were set aside.
#'
#' \strong{When to use.} Use \code{reliability_omega_categorical} when items are
#' ordered-categorical, especially when (a) the number of categories is
#' small (e.g., two to five), (b) item distributions are skewed, or (c)
#' threshold patterns differ markedly across items. In Kelley and
#' Pornprasertmanit's (2016) Study 3, treating ordered items as
#' continuous and using \code{\link{reliability_omega}} with the
#' observed total variance in the denominator achieved
#' acceptable coverage only when threshold patterns were similar across
#' items; in their experience that condition is rare in practice.
#'
#' Available confidence interval methods (set via \code{ci_method}).
#' Every interval here is bootstrap based: the rows of \code{data} are
#' resampled with replacement \code{B} times (10000 by default, 1000
#' under \code{missing = "mi"}; see \emph{Missing data} below) and
#' categorical omega is recomputed, with the full WLSMV model refit, on
#' each replication (Efron & Tibshirani, 1993). Replications whose
#' refit fails or does not converge, most common with small samples and
#' sparse response categories, are set aside, as are replications with
#' an improper solution, replications that leave a polychoric
#' correlation undetermined under the observed setting, and replications
#' whose determined correlations leave the loadings unidentified under
#' either setting (all described above), and the interval is
#' computed from the replications that return a value: the standard
#' error is their standard deviation, and the percentile and BCa limits
#' are quantiles of the returned replicates, not of all \code{B}. How
#' many were attempted, how many returned a value, and why the rest
#' were set aside is recorded in the \code{bootstrap} attribute of the
#' result (see \emph{Value}), printed in one line under the table, and
#' stated in a warning whenever a replication was set aside, so an
#' interval that conditions on a subset of the replications says so and
#' by how much. When fewer than half of the replications return a value,
#' no interval is reported: the standard error and the limits are
#' \code{NA}, and the warning and the printed line say why. The rest
#' would be a minority of the resamples, chosen by whether they returned
#' a value, and their spread would describe those resamples rather than
#' the sampling distribution of the coefficient; this is the criterion by
#' which the point estimate under \code{missing = "mi"} stops when fewer
#' than half of its completed data sets return a value (Lim & Cheung,
#' 2022, p. 1069). The remedy is a larger sample, which fills the sparse cells of
#' the bivariate tables. Bootstrap results vary from run to run; supply
#' \code{seed} for reproducibility.
#' \describe{
#'   \item{\code{"bca"}}{The bias-corrected and accelerated bootstrap,
#'   the default and the specific recommendation of Kelley and
#'   Pornprasertmanit (2016). Where the percentile interval reads its
#'   limits directly off the empirical quantiles of the bootstrap
#'   estimates, BCa adjusts the two quantile positions for median bias
#'   (estimated from the bootstrap distribution) and for the rate at
#'   which the estimator's variance changes with the parameter (the
#'   acceleration, estimated by the jackknife), making it second-order
#'   accurate where the percentile interval is first-order accurate
#'   (DiCiccio & Efron, 1996). The adjusted quantile positions sit
#'   farther into the tails than the percentile interval uses, which
#'   is why the default \code{B = 10000} is larger than the customary
#'   2000; reduce \code{B} for exploration, not for a reported
#'   analysis. The acceleration is estimated as \pkg{boot} estimates
#'   it: by regressing the bootstrap replicates on the resampling
#'   frequencies when at least as many replications as cases returned
#'   a value, at no extra cost, and otherwise, or when that regression
#'   fails, by the jackknife, one leave-one-out refit per case. A
#'   leave-one-out refit that fails (when a case is the only respondent
#'   in one of an item's categories, for example) is counted and set
#'   aside, and the
#'   acceleration is estimated from the refits that returned a value.
#'   When the acceleration cannot be estimated at all
#'   (fewer than two leave-one-out refits returned a value), or the
#'   bias correction cannot be (every returned replicate on one side of
#'   the estimate), the percentile interval is returned in its place
#'   with a warning, and the \code{ci_method} attribute records the
#'   substitution as \code{"percentile (BCa acceleration could not be
#'   estimated: ...)"}. The counts of leave-one-out refits and the
#'   acceleration used are in the \code{bootstrap} attribute.}
#'   \item{\code{"percentile"}}{Percentile bootstrap: the interval
#'   limits are the empirical quantiles of the bootstrap estimates.}
#'   \item{\code{"bootstrap_se"}, \code{"bootstrap_se_logistic"}}{Wald
#'   intervals using the bootstrap standard deviation as a standard
#'   error, built on the logit scale for the \code{_logistic} variant
#'   so the endpoints respect [0, 1]. Under \code{missing = "mi"} the
#'   standard error is that of the point estimate (see \emph{The number
#'   of imputations} below).}
#'   \item{\code{"none"}}{Return only the point estimate.}
#' }
#' Under \code{missing = "mi"} the percentile interval is the
#' recommended one and BCa is refused (see \emph{Missing data}).
#'
#' \strong{Missing data.} The \code{missing} argument sets how
#' incomplete rows of \code{data} are treated. Under \code{"pairwise"}
#' and \code{"mi"}, rows with no observed item, which carry no
#' information about the items, are set aside, as under full
#' information maximum likelihood elsewhere in the family, and when no
#' item value is missing in the rows that remain both return the
#' complete data result, with a message. Under multiple imputation a
#' row with no observed item is kept when an auxiliary variable is
#' observed on it, and its items are imputed from the auxiliaries, as
#' full information maximum likelihood elsewhere in the family keeps
#' such a row: setting it aside would make the analysis conditional on
#' having an observed item, which biases the estimate when the items go
#' missing together at random given the auxiliaries.
#' \describe{
#'   \item{\code{"listwise"} (default)}{Only the rows observed on every
#'   column of \code{data} are analyzed. Listwise deletion is consistent
#'   when values
#'   are missing completely at random (MCAR), and even then it discards
#'   what the incomplete rows say about the items they observe. When
#'   fewer than two rows are complete, or an item takes a single value in
#'   the complete rows, the function stops, says so, and names the
#'   treatment that uses the incomplete rows.}
#'   \item{\code{"pairwise"}}{Each item's thresholds are estimated from
#'   the rows that observe it and each polychoric correlation from the
#'   rows that observe both items (\pkg{lavaan}'s \code{missing =
#'   "pairwise"}); the single-factor fit, both denominators, and the
#'   assembly are otherwise unchanged. Pairwise deletion is the default
#'   treatment of missing data under WLSMV estimation in Mplus (Shi,
#'   Lee, Fairchild, & Maydeu-Olivares, 2020). Asparouhov and Muthén
#'   (2010) show it consistent when the missingness depends only on
#'   covariates in the model; a single-factor model of the items has no
#'   covariates, so here it is consistent only under MCAR, the
#'   assumption of listwise deletion. It uses more of the data than
#'   listwise deletion, not a weaker assumption: when values are missing
#'   at random (MAR), Shi et al. (2020) found the estimates of pairwise
#'   deletion severely biased, and they recommend multiple imputation.
#'   Every interval method applies, each replication resampling the rows
#'   and refitting by pairwise deletion.}
#'   \item{\code{"mi"}}{Multiple imputation, consistent when values are
#'   MAR given the items and any auxiliary variables named in
#'   \code{aux}. The items are imputed by chained equations with
#'   \pkg{mice} (van Buuren & Groothuis-Oudshoorn, 2011), with every
#'   other item and every auxiliary variable as a predictor and 10
#'   iterations of the chains, and categorical omega is computed on each
#'   completed data set with the chosen \code{denominator}; the point
#'   estimate is their mean (Rubin's rule for a point estimate), and the
#'   residual polychoric correlations are averaged the same way. An
#'   auxiliary variable is a predictor in the imputation model only,
#'   never an item, so it carries information about the missing values
#'   without entering the composite. The interval is the bootstrap then
#'   impute interval of Schomaker and Heumann (2018, their Boot MI):
#'   \code{B} bootstrap samples are drawn from the incomplete rows, each
#'   is imputed \code{imputations} times, the estimate of a replication
#'   is the mean of its imputations, and the percentile interval of
#'   those means is the confidence interval. Because the imputation
#'   model is refitted on every bootstrap sample, the interval carries
#'   the sampling variability of the imputation model as well as that of
#'   the coefficient. Bartlett and Hughes (2020) found that imputing and
#'   then bootstrapping does not, in general, give valid variances when
#'   the imputation and analysis models are uncongenial, as a regression
#'   imputation model and a probit factor model are, while bootstrapping
#'   and then imputing can. Schomaker and Heumann found the coverage of
#'   Boot MI near the nominal level in their simulations, and found that
#'   pooling the estimates of every imputation of every replication into
#'   one distribution gives intervals that are too wide, which
#'   \pkg{DMAR} therefore never does. A completed data set of the
#'   observed data that returns no value is set aside and counted, and
#'   the point estimate averages the rest; the call stops only when fewer
#'   than half return a value. A replication computes the same
#'   statistic: the mean of those of its completed data sets that return
#'   a value, set aside when fewer than half of them return one or when
#'   its imputation fails (see \emph{When a completed data set is
#'   improper} below).}
#' }
#'
#' \strong{When a completed data set is improper.} Imputation can produce
#' an improper single-factor solution (a standardized loading above one)
#' in a completed data set when the observed data have none, and the
#' example below shows how. The polychoric correlation of an item piled
#' into its top category with one piled into its bottom category is
#' determined by the few respondents in the cells where the two
#' disagree, and the values missing at random given \code{y1} and
#' \code{y2} are missing most often among the respondents low on them,
#' who supply those cells. Predictive mean matching fills each missing
#' value with the observed value of a donor whose predicted mean is
#' close, and among the observed respondents low on \code{y1} and
#' \code{y2} nearly all chose the bottom category of \code{y3} (all but
#' one, in the example's sample). How many imputed respondents land in
#' the rare cells therefore varies from one completed data set to the
#' next; where almost none do, the polychoric correlation runs toward
#' one and the loading of \code{y3} past it. In 100 completed data sets
#' of the example's data, 13 were improper: 7 of the 13 in which the
#' imputations put no respondent in the cells where \code{y3} is above
#' its bottom category and \code{y1} or \code{y2} below its top, 6 of
#' the 58 in which they put one or two there, and none of the 29 in which
#' they put three or more. The standard deviation of the loading of
#' \code{y3} over those imputations was 0.106, against 0.046 over 40
#' imputations under the proportional odds model. Neither that model
#' nor more donors is a remedy: over 61 samples of the example's design
#' the proportional odds model gave improper sets in fewer samples (5 of
#' the 61, against 13) but estimates lower than the complete data value
#' by 0.066 on average, and predictive mean matching with ten donors in
#' place of five lowered the average share of improper sets only from
#' 3\% to 2\%.
#'
#' Categorical omega is not defined for an improper solution, so such a
#' completed data set has no estimate, and Rubin's rules, which average
#' an estimate from every completed data set, need one. Stopping
#' whenever one completed data set fails, as Jia and Wu (2019) counted a
#' replication of their simulation as failed unless all its imputations
#' converged, makes the method fail by design: with a share \eqn{p} of
#' such sets, \eqn{M} imputations include one with probability
#' \eqn{1 - (1 - p)^M}, which approaches one as von Hippel's rule asks
#' for more imputations. The function instead sets such a completed
#' data set aside, as it sets aside a bootstrap replication, and the
#' point estimate averages those that returned a value, as Lim and
#' Cheung (2022) pooled the imputations that converged and as
#' \pkg{lavaan.mi} does by default for imputations that did not
#' converge. Setting aside conditions the estimate on a property of the
#' estimate, so the number set aside, and why, is reported in a warning
#' and printed under the table, and the call stops when fewer than half
#' of the completed data sets return a value, the criterion by which Lim
#' and Cheung counted an analysis as failed. The improper sets are those
#' with the largest loadings, so setting them aside lowers the estimate.
#' Over the 61 samples, 13 had at least one improper completed data set
#' among 20, with a share of 3\% on average and 65\% at most, and no
#' sample's complete data or pairwise deletion fit was improper. Against
#' refitting those sets with every loading bounded at one, setting them
#' aside lowered the estimate by about 0.034 times the share set aside
#' (0.0017 at 5\%, a tenth of the complete data estimate's sampling
#' standard deviation of 0.0165), and over the 13 samples the two were
#' equally far from the complete data estimate (root mean square
#' differences 0.0126 and 0.0122). No completed data set of the
#' observed data is fitted after the bootstrap, so none can discard a
#' bootstrap already run (see \emph{The number of imputations}).
#'
#' A bootstrap replication follows the same rule, set by set: it
#' averages those of its \eqn{M_b} completed data sets that return a
#' value and is set aside when fewer than half do, because the bootstrap
#' distribution is that of the statistic computed as the point estimate
#' is computed. The interval follows it too: when fewer than half of the
#' replications return a value, none is reported, as for the other
#' treatments (see the confidence interval methods above). Setting a replication aside whenever any one of its
#' completed data sets returned no value would be a different statistic,
#' and the argument above applies to it unchanged: a replication would
#' be set aside with probability \eqn{1 - (1 - p)^{M_b}}, so raising
#' \code{imputations} would set more replications aside, and the
#' completed data sets of those replications that did return a value
#' would be discarded. In a bootstrap of the example's incomplete data
#' under the model implied setting (\code{B = 200},
#' \code{imputations = 5}, \code{seed = 1}), in which a quarter of the
#' completed data sets returned no value, that rule, applied to the first
#' two, three, four, and five imputations of each replication, would have
#' set aside 35.5, 45.5, 48.5, and 52.5 percent of the replications,
#' where the rule used sets aside 12.5, 21, 14.5, and 20 percent. Under
#' the rule used the share set aside goes to zero as \eqn{M_b} grows
#' when fewer than half of the completed data sets fail, but not
#' steadily: an odd \eqn{M_b} needs \eqn{(M_b + 1)/2} of its sets to
#' return a value, one more than the even number below it needs of one
#' set fewer, so it sets aside more replications than either even
#' neighbor (with a fifth of the sets failing independently, 4.0, 10.4,
#' 2.7, and 5.8 percent at two to five imputations). Where completed
#' data sets are set aside, an even number of imputations is the better
#' choice, and the number the warning names is even. The warning and the
#' \code{imputation} attribute count
#' the completed data sets set aside within the replications kept
#' (\code{replicate_sets_set_aside}).
#'
#' \strong{The imputation model.} \code{impute = "pmm"}, the default, is
#' predictive mean matching: each missing value is replaced by the
#' observed value of the same item of a donor whose predicted mean, from
#' a linear regression on the other variables, is close, with the
#' regression parameters drawn anew for each imputation. It never
#' imputes a category the item does not show, though a category rare
#' among the observed values is rarely drawn. \code{impute = "polr"} is
#' the proportional odds model, with every other item entering as a
#' factor (a binary item is imputed by logistic regression, the same
#' model for two categories). \pkg{mice} draws its categories at the
#' model's point estimates, without drawing the model's parameters, and
#' Jia and Wu (2019), pooling its imputations in the usual way, found
#' acceptable loading estimates but confidence intervals that were too
#' narrow under missingness at random; under bootstrap then impute the
#' refit on every bootstrap sample supplies the variability of the
#' imputation model instead. The default rests on a Monte Carlo simulation
#' run for this function: 10 items, one factor with loadings from .50
#' to .80, \eqn{n = 300}, and items 6 to 10 missing at random given the
#' mean of items 1 to 5, about 21\% of the values in all. Against the
#' complete data estimate of the same sample, the mean categorical omega
#' was off by -0.0035 with five categories (100 replications) and by
#' +0.0015 with three skewed categories (60 replications) under
#' predictive mean matching; by -0.019 and +0.035 under the
#' proportional odds model; by -0.028 and -0.047 under pairwise
#' deletion; and by -0.105 with five categories under listwise deletion,
#' which with three categories failed in 55 of the 60 samples.
#' A replication under predictive mean matching was also about ten
#' times faster (see \emph{Run time} below). The
#' evidence is one factor structure and one missing data mechanism, so
#' with very rare categories or another mechanism the comparison could
#' differ.
#'
#' Every incomplete item is imputed. By default \pkg{mice} stops
#' imputing a variable whose correlation with another, computed from
#' the rows that observe both, reaches .999 in absolute value; two items
#' observed together on only two or three rows, as in a planned missing
#' design with a small overlap, reach it by arithmetic, and the item
#' left with its missing values would then be fitted by listwise
#' deletion inside each completed data set. That check is therefore
#' turned off, an item with a single observed value is refused before
#' anything is imputed, and an item \pkg{mice} still does not impute
#' stops the call with its name.
#'
#' The same check also set aside an auxiliary variable that duplicates
#' another (a score and its rescaling, age in years and in months).
#' Without it both reach every item's regression, and when the
#' predictors carry a second exact dependency (a second such pair, or an
#' auxiliary variable computed from complete items) \pkg{mice} 3.19.0
#' removes only one of them and falls back to a ridge penalty, after
#' which its imputations are unstable: on the five-item example of the
#' tests, an auxiliary variable computed from the first three items and
#' its rescaling gave 0.768 with a between-imputation standard deviation
#' of 0.062, against 0.800 and 0.008 without them. So the auxiliary
#' variables are screened before anything is imputed: one whose multiple
#' correlation with the columns observed on every row that observes it
#' (the items and the auxiliary variables kept before it) reaches .999,
#' \pkg{mice}'s own threshold for a pair, is a linear function of
#' predictors already in the model, carries no information about the
#' missing values, and is left out, with a message naming it; the items
#' are not screened, so a planned missing design is still imputed in
#' full. When \pkg{mice} still applies its ridge penalty, which it does
#' when the predictors of an item are (nearly) linearly dependent, the
#' warning says how often, and the \code{imputation} attribute counts it
#' (\code{ridge}).
#'
#' The items are imputed as plain integer codes, whatever their class: a
#' column read with \pkg{haven}, which carries value labels, and a
#' factor whose labels are whole numbers are both accepted (see the
#' \code{data} argument), and the imputation model is the one for the
#' codes. SPSS user-defined missing values, which
#' \code{haven::read_sav(user_na = TRUE)} keeps as values of a
#' \code{haven_labelled_spss} column, are missing under every treatment,
#' as SPSS and \pkg{haven}'s \code{is.na()} treat them.
#'
#' \strong{The number of imputations.} Each bootstrap replication is
#' imputed \code{imputations} times (\eqn{M_b}, 20 by default). The
#' variance of a replicate mean is the sampling variance of the
#' coefficient plus the imputation variance divided by \eqn{M_b}, so
#' imputation noise widens the percentile interval by about
#' \eqn{\sqrt{1 + \gamma / M_b}}, where \eqn{\gamma} is the fraction of
#' missing information: by 1.5\% at \eqn{\gamma = .15} and by 4.9\% at
#' \eqn{\gamma = .50} with five imputations. The function estimates
#' \eqn{\gamma} from the replications themselves, by the one-way
#' analysis of variance of the replicate by imputation estimates (the
#' variance between replications net of the imputation variance within
#' them, which is the sampling variance with infinitely many
#' imputations, against the imputation variance), reports the widening,
#' and warns when it exceeds 5\%, naming the value of \code{imputations}
#' that would bring it under; when no interval is reported (fewer than
#' half of the replications returned a value) neither the warning nor the
#' printed line gives a widening, and the printed fraction of missing
#' information is marked as that of the replications that returned one. A replication that sets some of its
#' completed data sets aside averages \eqn{n_b < M_b} of them, and the
#' analysis of variance is then the unbalanced one (Searle, Casella, &
#' McCulloch, 1992, Sec. 3.6), whose interval for the ratio of the
#' variances is approximate; the widening is
#' \eqn{\sqrt{1 + \gamma \, \overline{1/n_b}}}, the mean taken over the
#' replications, and the value of \code{imputations} named allows for
#' the share of completed data sets set aside and is even when any were,
#' whether in the replications kept or in those set aside for them, but
#' not when only the imputation itself failed (see \emph{When a completed
#' data set is improper}). Schomaker and Heumann (2018) found that
#' Boot MI can do well with fewer than five imputations per replication
#' and covers too often with one. The point estimate averages
#' \eqn{M_0} imputations of the observed data: at least 20 and at least
#' \code{imputations}, and more when von Hippel's (2020, Eq. 1)
#' quadratic rule, \eqn{M = 1 + \frac{1}{2}(\gamma / \mathrm{CV})^2} with
#' \eqn{\mathrm{CV} = .05}, asks for more at the upper limit of a 95\%
#' confidence interval for \eqn{\gamma}, as his two stage procedure
#' uses it (here the limit for the ratio of variance components in the
#' same analysis of variance). His target is the replicability of a
#' standard error from Rubin's rules, which bootstrap then impute does
#' not use; the rule is borrowed to keep the point estimate from
#' changing much if the data were imputed again, with the bootstrap as
#' the pilot. Since the rule asks for at most 201 imputations (at
#' \eqn{\gamma = 1}), a call with a bootstrap makes that many
#' imputations of the observed data, and fits and checks each, before
#' the bootstrap runs, and the point estimate then uses the first
#' \eqn{M_0} of them, the number fixed by the rule and the floor and
#' never by their estimates. Nothing found in a completed data set of
#' the observed data can then discard a finished bootstrap; the
#' imputations the rule does not use cost at most 181 imputations and
#' fits, less than 1\% of the 20000 of the default bootstrap. The estimate of
#' \eqn{\gamma} is imprecise
#' when \code{B} is small, and under \code{impute = "polr"}, whose
#' imputations ignore the uncertainty of the imputation model's
#' parameters, it is too low. Without a bootstrap the point estimate
#' averages \code{max(imputations, 20)} imputations, which keeps Rubin's
#' relative efficiency, \eqn{(1 + \gamma / M)^{-1}}, at .976 or above for
#' \eqn{\gamma} up to .5. The \code{se} row is the standard error of the
#' point estimate, \eqn{\sqrt{s^2_{boot} + s^2_{imp} / M_0}}, with
#' \eqn{s^2_{boot}} the sampling variance and \eqn{s^2_{imp}} the
#' imputation variance from the same analysis of variance, so it counts
#' the noise of the point estimate's own imputations. The analysis of
#' variance is a construction for this function, not a published
#' procedure for omega.
#'
#' \strong{Run time under multiple imputation.} A replication imputes
#' its resample \code{imputations} times and refits the factor model on
#' each completed data set, so the run time grows with \code{B} times
#' \code{imputations}. With 10 items, five categories, and 300 cases, a
#' replication with five imputations took about half a second under
#' predictive mean matching and about five seconds under the
#' proportional odds model, so the default \code{B = 1000} runs for
#' about eight minutes and well over an hour, respectively, on one
#' core with five imputations. The default of 20 imputations per
#' replication does about four times that imputation and fitting work.
#' This larger default reduces imputation noise; it has not been
#' separately calibrated for interval coverage. The default \code{B} is 1000 under
#' \code{missing = "mi"} rather than the 10000 of the other two
#' treatments (a percentile interval reads its limits from less extreme
#' quantiles than BCa), and why BCa, whose jackknife would add
#' \code{imputations} imputations per case, is refused. \code{ncpus}
#' spreads the replications over processes; every resample, and a seed
#' for each replication's imputations, is drawn before any replication
#' runs, so a given \code{seed} returns the same result whatever
#' \code{ncpus} is.
#'
#' \strong{Comparison with other packages.} The \pkg{psych} package's
#' \code{\link[psych]{omega}} fits a Schmid-Leiman hierarchical factor
#' model to Pearson correlations or, with \code{poly = TRUE}, to
#' polychoric correlations, and does not implement categorical omega in
#' the sense of Green and Yang (2009), which
#' \code{reliability_omega_categorical} computes for ordered-categorical
#' items; \code{\link[psych]{polychoric}} provides polychoric correlation
#' estimation as a separate tool.
#'
#' @param data A numeric matrix or data frame of three or more
#'   ordered-categorical items (integer codes for the categories), plus
#'   any auxiliary variables named in \code{aux}. Only the order of an
#'   item's codes matters: its categories are scored 0, 1, \ldots,
#'   \eqn{C - 1} (see \emph{How the categories are scored} in
#'   \emph{Details}). A factor whose labels are whole numbers is read as
#'   those numbers, so a factor labeled 0 to 4 is the item coded 0 to 4,
#'   under every treatment of \code{missing}; a factor whose labels are
#'   words is refused, as is an ordered factor whose levels are not in the
#'   order of the numbers they name. A column read with \pkg{haven} is
#'   read as its values, with SPSS user-defined missing values made
#'   missing. Incomplete rows are
#'   treated as \code{missing} says; under the default they are set
#'   aside. With two items the single-factor model is not identified
#'   (their one polychoric correlation determines only the product of
#'   the two loadings), and the function stops.
#' @param ci_method Method for constructing the confidence interval. See
#'   \emph{Details}. When not supplied, no interval is computed: every
#'   interval for categorical omega is bootstrap based, and a bootstrap
#'   is never run unless requested. Ask for \code{"bca"} (the
#'   recommended method) or \code{"percentile"} to obtain one.
#' @param conf_level Confidence level for the interval. Defaults to
#'   \code{0.95}.
#' @param B Number of bootstrap replications. The default, \code{NULL},
#'   is 10000 under \code{missing = "listwise"} and \code{"pairwise"},
#'   where BCa, the recommended interval, reads its limits from extreme
#'   quantiles, and 1000 under \code{missing = "mi"}, where each
#'   replication is imputed \code{imputations} times and refitted on
#'   each, and the percentile interval is the recommended one (see
#'   \emph{Run time under multiple imputation}). When no item value is
#'   missing, \code{missing = "mi"} has nothing to impute and returns
#'   the complete data result, whose bootstrap refits once per
#'   replication, so the default there is 10000, as the message says.
#'   When a bootstrap is requested, \code{B} (supplied or
#'   default) must be a single whole number of at least the count the
#'   interval method and the level require, as in the rest of the family:
#'   at the default \code{conf_level}, 400 for \code{"percentile"} (10
#'   replications beyond each limit), 2000 for \code{"bca"} (50 beyond
#'   each nominal limit), and 100 for the two bootstrap standard error
#'   methods; at 0.99 the first two are 2000 and 10000, so the default of
#'   1000 under \code{missing = "mi"} does not serve a 99 percent
#'   percentile interval and the refusal names the number. A smaller
#'   \code{B} is refused before anything is computed (class
#'   \code{dmar_count_invalid}). When fewer replications than required return
#'   a value the interval is withheld under a warning of class
#'   \code{dmar_bootstrap_insufficient}. See
#'   \code{\link{dmar_bootstrap}}, which also describes
#'   \code{options(dmar.bootstrap_exploratory = TRUE)} for a first look
#'   with a smaller \code{B}.
#' @param seed Random number seed used for bootstrap reproducibility.
#'   Defaults to \code{NULL}, under which the resamples are drawn from,
#'   and advance, the current random number stream; supply an integer
#'   for a reproducible interval, which is set for the call only, the
#'   caller's random number state being restored when it returns.
#' @param denominator How the polychoric correlations behind the total
#'   variance in the denominator are estimated: \code{"observed"} (the
#'   default; the unrestricted polychoric correlations) or
#'   \code{"model_implied"} (the correlations reproduced by the fitted
#'   single-factor model, as in Green and Yang, 2009). See
#'   \emph{Details}. Every interval method applies to both settings,
#'   since every interval here is bootstrap based and each replication
#'   recomputes the coefficient with the same denominator.
#' @inheritParams reliability_omega
#' @param ncpus The number of processes over which the bootstrap
#'   replications are spread. The default \code{1} runs them in this R
#'   session; a larger number uses forked workers (\pkg{boot}'s
#'   \code{parallel = "multicore"}) or, on Windows, which cannot fork, a
#'   socket cluster (\code{"snow"}), whose workers load the installed
#'   \pkg{DMAR}. \pkg{boot} draws every resample before any replication
#'   runs, so a given \code{seed} returns the same interval, counts, and
#'   replicates whatever \code{ncpus} is. A session that forks workers and also
#'   starts background processes through \pkg{callr} or \pkg{processx}
#'   needs the environment variable \code{PROCESSX_NOTIFY_OLD_SIGCHLD} set
#'   before \pkg{processx} is loaded (see \code{\link{sem_exact_fit_bootstrap}}). Ignored unless a bootstrap
#'   method is selected.
#' @param missing How incomplete rows of \code{data} are treated:
#'   \code{"listwise"} (the default; only complete rows are analyzed),
#'   \code{"pairwise"} (each threshold and
#'   polychoric correlation from the rows that observe its items), or
#'   \code{"mi"} (multiple imputation, with a bootstrap then impute
#'   interval). Listwise and pairwise deletion are consistent only when
#'   values are missing completely at random; multiple imputation when
#'   they are missing at random given the items and \code{aux}. See
#'   \emph{Missing data} in \emph{Details}.
#' @param aux Optional character vector naming columns of \code{data}
#'   that are auxiliary variables: numeric predictors in the imputation
#'   model, never items, so they inform the imputations without entering
#'   the composite (the inclusive strategy of Collins, Schafer, & Kam,
#'   2001). Supplying \code{aux} implies \code{missing = "mi"}, as it
#'   implies full information maximum likelihood in
#'   \code{\link{reliability_omega}}; with \code{"listwise"} or
#'   \code{"pairwise"} named explicitly it is an error.
#' @param impute The imputation model under \code{missing = "mi"}:
#'   \code{"pmm"} (the default; predictive mean matching) or
#'   \code{"polr"} (the proportional odds model). See \emph{The
#'   imputation model} in \emph{Details}. Supplying it with another
#'   \code{missing} is an error.
#' @param imputations Under \code{missing = "mi"}, the number of
#'   imputations of each bootstrap replication (at least 2, default 20),
#'   which is also the least number of imputations behind the point
#'   estimate (whose number is at least 20 and may be raised by von
#'   Hippel's rule; see \emph{The number of imputations}). Supplying it
#'   with another \code{missing} is an error.
#'
#' @return A \code{data.frame} with columns \code{term} and \code{value}
#' and rows
#' \code{"estimate"} (sample \eqn{\omega_C}),
#' \code{"se"} (the bootstrap standard deviation across replications,
#' already on the coefficient scale; \code{NA} for
#' \code{ci_method = "none"}),
#' \code{"ci_lower"} and \code{"ci_upper"} (clamped to [0, 1]),
#' \code{"conf_level"}, \code{"N"} (the rows analyzed: the complete
#' rows under listwise deletion, otherwise the rows with at least one
#' observed item, and under multiple imputation with \code{aux} also the
#' rows with an observed auxiliary variable), \code{"N_complete"} (the
#' rows complete on the items;
#' equal to \code{"N"} under listwise deletion, and carried so the whole
#' reliability family returns one shape), and \code{"n_items"} (number of items); the rows are
#' the same under every \code{missing}. Under \code{missing = "mi"} the
#' estimate is the mean over the imputations of the observed data and
#' \code{se} is the standard error of that mean from the bootstrap (see
#' \emph{The number of imputations}). Attributes \code{coefficient}
#' (\code{"omega_categorical"}), \code{denominator} (\code{"observed"}
#' or \code{"model_implied"}), \code{missing} (\code{"listwise"},
#' \code{"pairwise"}, or \code{"mi"}; \code{"listwise"} also when the
#' items are complete and another treatment had nothing to do),
#' \code{aux} (when supplied), \code{ci_method} (under
#' \code{missing = "mi"} the method followed by \code{"(bootstrap then
#' impute)"}), and \code{B} record the computation. The attribute
#' \code{residual_cor} carries the local dependence summary of the
#' single-factor fit: a list
#' with \code{matrix} (the unrestricted polychoric correlations minus
#' those the model reproduces, named by item), \code{largest} (the
#' residual correlation largest in absolute value, with its sign),
#' \code{items} (the pair it belongs to), and \code{saturated}; a pair
#' whose polychoric correlation the data do not determine has an
#' \code{NA} residual and is passed over. With three items the
#' single-factor model is saturated (zero degrees of freedom) and
#' reproduces every polychoric correlation, so its residuals are zero
#' and carry no information about local dependence: \code{saturated} is
#' then \code{TRUE}, \code{largest} is \code{NA}, and \code{items} is
#' \code{NULL}, rather than naming whichever pair carries the largest
#' rounding error. The printed table names the denominator and the
#' largest residual correlation (or says that the model is saturated)
#' in two lines under the table. A bootstrap call also carries the attribute
#' \code{bootstrap}, a list recording the replications the interval
#' conditions on: \code{attempted} (\code{B}), \code{returned} (the
#' replications whose value the interval uses; the standard error is
#' their standard deviation and the percentile and BCa limits are their
#' quantiles), \code{failed} (the difference), \code{reasons} (a named
#' count of why replications were set aside: \code{error} when the refit
#' raised an error, \code{non_finite} when it was set aside without one,
#' as a fit that \pkg{lavaan} stopped or that did not converge is, as an
#' improper one is, and, under the observed denominator, as one that
#' left a polychoric correlation undetermined is, and, under either, as
#' one whose determined correlations leave the loadings unidentified is),
#' \code{error_messages} (the
#' distinct error messages with their counts, empty when none),
#' \code{undetermined} (the replications set aside because the resample
#' left a polychoric correlation undetermined, which the observed
#' setting needs, or left correlations undetermined whose absence leaves
#' the loadings unidentified, under either setting; see
#' \emph{Details}), \code{improper} (the replications
#' set aside because their single-factor solution was improper, with a
#' standardized loading above one, which is an error variance below
#' zero), \code{improper_kept} (always zero here: categorical omega is
#' not defined for an improper solution, so such a replication has no
#' estimate to keep), \code{improper_rule} (\code{"drop"}, for the same
#' reason), \code{no_solution} (always zero here: it counts the
#' replications that \code{\link{reliability_omega}} sets aside because
#' the correlations of three items admit no single-factor solution, a
#' check the categorical fit does not make),
#' \code{not_converged} (zero here: a fit that did not converge is
#' counted among the \code{non_finite} ones),
#' \code{too_few_cases} (always zero here: it counts the replications
#' that \code{\link{reliability_omega}} sets aside for too few distinct
#' cases, a rule the categorical fit does not apply),
#' \code{caution} (\code{TRUE} when the replications set aside number
#' at least \eqn{(1 - C)/2} of \code{B}, with \eqn{C} the
#' confidence level, so that where they would have fallen can decide a
#' limit), \code{require_half} (\code{TRUE}: no interval is formed when
#' fewer than half of the replications return a value),
#' \code{required}, \code{bootstrap_status}, and \code{B_used}
#' (the fewest replications a request may ask for, the verdict of the rule,
#' and the count the limits rest on, as in \code{\link{reliability_alpha}}),
#' \code{acceleration} and \code{acceleration_method} (the BCa
#' acceleration used and whether it came from the regression or the
#' jackknife; \code{NA} unless a BCa interval was formed),
#' \code{jackknife} (the same \code{attempted}, \code{returned},
#' \code{failed}, \code{reasons}, \code{error_messages},
#' \code{undetermined}, \code{improper}, \code{improper_kept},
#' \code{no_solution}, \code{not_converged}, and
#' \code{too_few_cases} for the leave-one-out refits behind the
#' acceleration, all zero when no jackknife ran), \code{ncpus} (the
#' number of processes used), \code{t} (the estimate of every
#' replication, in order, \code{NA}
#' where there is none), and \code{status} (a factor giving the fate of
#' every replication, with levels \code{"proper"}, \code{"improper"},
#' \code{"error"}, \code{"non_finite"}, \code{"undetermined"}, and
#' \code{"not_converged"}; the interval uses those whose status is
#' \code{"proper"}). The attribute
#' \code{B_used} is the number of replications the interval uses
#' (\code{returned}; when fewer than half returned a value, or fewer
#' than \code{required}, no interval is formed). The attributes
#' \code{n_dropped}, \code{B_required}, and \code{bootstrap_status}
#' repeat the counts in the names every bootstrap result of the package
#' carries (see \code{\link{dmar_bootstrap}}). The printed table shows these counts in one line,
#' and a warning states them whenever a replication was set aside or
#' the caution applies.
#'
#' Under \code{missing = "mi"} the \code{bootstrap} attribute has the
#' same elements, with three differences: \code{reasons} has a third
#' count, \code{imputation}, for the replications whose imputation
#' itself failed (not part of \code{non_finite}), whose messages are in
#' \code{imputation_messages}; \code{status} has a seventh level,
#' \code{"imputation"}; and a replication counted under another reason
#' is one in which fewer than half of its completed data sets returned a
#' value, the first that returned none deciding the reason. The jackknife
#' counts are zero, since BCa is refused. The attribute
#' \code{imputation} is a list: \code{impute};
#' \code{imputations_point} (\eqn{M_0}, the imputations of the observed
#' data behind the point estimate); \code{imputations_returned} (those
#' of them that returned a value, which the estimate averages);
#' \code{imputations_set_aside} (the others, counted by reason:
#' \code{improper}, \code{undetermined}, \code{error},
#' \code{non_converged}, and \code{non_finite});
#' \code{imputations_replicate} (\code{imputations});
#' \code{imputations_floor}
#' (\code{max(imputations, 20)}); \code{imputations_rule} (the number
#' von Hippel's rule asks for); \code{maxit} (10, the iterations of the
#' chained equations); \code{estimates} (the \eqn{M_0} completed data
#' estimates, \code{NA} where one was set aside), \code{status} (the
#' fate of each: \code{"proper"} or the reason it was set aside), and
#' \code{between_sd} (the standard deviation of those that returned a
#' value);
#' \code{sd_replicates} (the standard deviation of the replicate
#' means); \code{s2_boot}, \code{s2_imp}, \code{fmi}, and
#' \code{fmi_ci_upper} (the sampling and imputation variances, the
#' fraction of missing information, and its 95\% upper limit);
#' \code{widening} (\eqn{\sqrt{1 + \gamma / M_b}}, the factor by which
#' imputation noise widens the percentile interval);
#' \code{replicate_estimates} (the \code{B} by \code{imputations}
#' matrix of completed data estimates, \code{NA} where one returned no
#' value; \code{NULL} without a bootstrap);
#' \code{replicate_sets_set_aside} (the completed data sets set aside
#' within the replications that returned a value, counted by the reasons
#' of \code{imputations_set_aside}; \code{NULL} without a bootstrap);
#' \code{polr_fallbacks} (how often the proportional odds model fell
#' back to a multinomial model, in the point estimate's imputations and
#' in the replications); \code{ridge} (how often \pkg{mice} applied its
#' ridge penalty, counted the same way); and \code{aux_set_aside} (the
#' auxiliary variables left out of the imputation model as linear
#' functions of other predictors, \code{character(0)} when none).
#' The components that need the bootstrap are \code{NA} without one.
#' Under \code{missing = "pairwise"} or \code{"mi"} the printed table
#' adds a line naming the treatment and its assumption, and under
#' \code{"mi"} a line giving the number of imputations, the fraction of
#' missing information, and the widening (none when no interval is
#' reported).
#'
#' @references
#' Asparouhov, T., & Muthén, B. (2010). \emph{Weighted least squares
#'   estimation with missing data} (Mplus technical appendix). Muthén &
#'   Muthén.
#'
#' Bartlett, J. W., & Hughes, R. A. (2020). Bootstrap inference for
#'   multiple imputation under uncongeniality and misspecification.
#'   \emph{Statistical Methods in Medical Research, 29}(12), 3533--3546.
#'   \doi{10.1177/0962280220932189}
#'
#' Collins, L. M., Schafer, J. L., & Kam, C.-M. (2001). A comparison of
#'   inclusive and restrictive strategies in modern missing data
#'   procedures. \emph{Psychological Methods, 6}(4), 330--351.
#'   \doi{10.1037/1082-989X.6.4.330}
#'
#' DiCiccio, T. J., & Efron, B. (1996). Bootstrap confidence intervals.
#'   \emph{Statistical Science, 11}(3), 189--228.
#'
#' Drezner, Z., & Wesolowsky, G. O. (1990). On the computation of the
#'   bivariate normal integral. \emph{Journal of Statistical Computation
#'   and Simulation, 35}(1--2), 101--107.
#'
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Genz, A. (2004). Numerical computation of rectangular bivariate and
#'   trivariate normal and t probabilities. \emph{Statistics and
#'   Computing, 14}(3), 251--260.
#'   \doi{10.1023/B:STCO.0000035304.20635.31}
#'
#' Green, S. B., & Yang, Y. (2009). Reliability of summed item scores
#'   using structural equation modeling: An alternative to coefficient
#'   alpha. \emph{Psychometrika, 74}, 155--167.
#'   \doi{10.1007/s11336-008-9099-3}
#'
#' Jia, F., & Wu, W. (2019). Evaluating methods for handling missing
#'   ordinal data in structural equation modeling. \emph{Behavior
#'   Research Methods, 51}, 2337--2355. \doi{10.3758/s13428-018-1187-4}
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
#' Lim, A. J.-M., & Cheung, M. W.-L. (2022). Evaluating FIML and
#'   multiple imputation in joint ordinal-continuous measurements models
#'   with missing data. \emph{Behavior Research Methods, 54}, 1063--1077.
#'   \doi{10.3758/s13428-021-01582-w}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' Millsap, R. E., & Yun-Tein, J. (2004). Assessing factorial invariance
#'   in ordered-categorical measures.
#'   \emph{Multivariate Behavioral Research, 39}(3), 479--515.
#'   \doi{10.1207/s15327906mbr3903_4}
#'
#' Muthén, B. (1984). A general structural equation model with
#'   dichotomous, ordered categorical, and continuous latent variable
#'   indicators. \emph{Psychometrika, 49}(1), 115--132.
#'
#' Schomaker, M., & Heumann, C. (2018). Bootstrap inference when using
#'   multiple imputation. \emph{Statistics in Medicine, 37}(14),
#'   2252--2266. \doi{10.1002/sim.7654}
#'
#' Searle, S. R., Casella, G., & McCulloch, C. E. (1992). \emph{Variance
#'   components}. Wiley.
#'
#' Shi, D., Lee, T., Fairchild, A. J., & Maydeu-Olivares, A. (2020).
#'   Fitting ordinal factor analysis models with missing data: A
#'   comparison between pairwise deletion and multiple imputation.
#'   \emph{Educational and Psychological Measurement, 80}(1), 41--66.
#'   \doi{10.1177/0013164419845039}
#'
#' Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
#'   reliability coefficients: Accuracy in parameter estimation via narrow
#'   confidence intervals. \emph{British Journal of Mathematical and
#'   Statistical Psychology, 65}, 371--401.
#'   \doi{10.1111/j.2044-8317.2011.02030.x}
#'
#' van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate
#'   imputation by chained equations in R. \emph{Journal of Statistical
#'   Software, 45}(3), 1--67. \doi{10.18637/jss.v045.i03}
#'
#' von Hippel, P. T. (2020). How many imputations do you need? A
#'   two-stage calculation using a quadratic rule. \emph{Sociological
#'   Methods & Research, 49}(3), 699--718.
#'   \doi{10.1177/0049124117747303}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{reliability}} (general wrapper),
#' \code{\link{reliability_omega}} (use for continuous items),
#' \code{\link{reliability_kr20}} (for dichotomous items),
#' \code{\link[psych]{omega}},
#' \code{\link[psych]{polychoric}}.
#'
#' @examples
#' set.seed(113)
#' # The mixed-skew design of Green and Yang (2009, Table 1): four items
#' # loading 0.8 on one factor, each cut into five categories, two piled
#' # into the top category and two into the bottom. Their population
#' # categorical omega for this design is 0.716.
#' N <- 1000
#' p_top <- c(.011, .027, .047, .118, .797)
#' cut_top <- c(-Inf, qnorm(cumsum(p_top)[-5]), Inf)
#' cut_bottom <- c(-Inf, qnorm(cumsum(rev(p_top))[-5]), Inf)
#' latent <- outer(rnorm(N), rep(0.8, 4)) +
#'           matrix(rnorm(N * 4, sd = 0.6), N, 4)
#' items <- cbind(y1 = cut(latent[, 1], cut_top, labels = FALSE),
#'                y2 = cut(latent[, 2], cut_top, labels = FALSE),
#'                y3 = cut(latent[, 3], cut_bottom, labels = FALSE),
#'                y4 = cut(latent[, 4], cut_bottom, labels = FALSE))
#'
#' # Categorical omega as Green and Yang defined it: the total variance
#' # of the summed categories comes from the correlations the fitted
#' # single-factor model reproduces.
#' reliability_omega_categorical(data = items, denominator = "model_implied")
#'
#' # The default, categorical omega as Kelley and Pornprasertmanit (2016)
#' # computed and recommended it: the total variance comes from the
#' # unrestricted polychoric correlations. The single factor holds in
#' # this design, so the two agree to within about 0.01; they part ways
#' # as the model misfits.
#' reliability_omega_categorical(data = items)
#'
#' # The same sum score with the items treated as continuous. Robust
#' # omega falls far below 0.716 here, because items whose categories
#' # pile up at opposite ends relate to the factor far from linearly,
#' # which is the case Green and Yang designed their coefficient for.
#' reliability_omega(data = items)
#'
#' # The same items with values missing at random: respondents low on y1
#' # and y2 leave y3 and y4 unanswered more often. Listwise deletion, the
#' # default, keeps only the complete rows, and here that leaves too
#' # little: the few respondents who determine the correlation of a
#' # top-piled item with a bottom-piled one (one to five of the 1000 for
#' # each pair) stand low on y1 and y2, the complete rows keep none of
#' # them for y2 with y3 or for y1 with y4, and the correlations left do
#' # not identify the loadings, so a call with the default stops and says
#' # so (see Details).
#' z <- as.numeric(scale(items[, "y1"] + items[, "y2"]))
#' items_mar <- items
#' items_mar[runif(N) < plogis(-1.5 - 1.2 * z), "y3"] <- NA
#' items_mar[runif(N) < plogis(-1.5 - 1.2 * z), "y4"] <- NA
#'
#' # Pairwise deletion estimates each polychoric correlation from the rows
#' # that observe both items, so every row contributes. It still assumes
#' # values missing completely at random, as listwise deletion does.
#' # Multiple imputation, which assumes only that they are missing at
#' # random, takes longer than an example here can (see Details).
#' reliability_omega_categorical(data = items_mar,
#'                               denominator = "model_implied",
#'                               missing = "pairwise")
#'
#' @keywords htest multivariate
#' @family reliability
#'
#' @export

reliability_omega_categorical <- function(data,
                                ci_method = c("bca", "percentile",
                                              "bootstrap_se",
                                              "bootstrap_se_logistic",
                                              "none"),
                                conf_level = 0.95,
                                B = NULL,
                                seed = NULL,
                                denominator = c("observed",
                                                "model_implied"),
                                ncpus = 1L,
                                missing = c("listwise", "pairwise", "mi"),
                                aux = NULL,
                                impute = c("pmm", "polr"),
                                imputations = 20L,
                                progress = NULL) {
  .dmar_check_conf_level(conf_level)
  ci_supplied <- !missing(ci_method)
  ci_method <- match.arg(ci_method)
  denominator <- match.arg(denominator)
  ncpus <- .check_ncpus(ncpus)
  progress <- .check_progress(progress)
  missing_supplied <- !missing(missing)
  missing <- match.arg(missing)
  impute_supplied <- !missing(impute)
  imputations_supplied <- !missing(imputations)
  impute <- match.arg(impute)

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  }
  if (missing(data) || is.null(data)) {
    .dmar_refuse_argument(c(
      "'data' (a matrix or data frame of ordered-categorical item ",
      "scores) is required."), arg = "data")
  }
  # When a bootstrap is requested, B follows the rule of ?dmar_bootstrap
  # (.check_boot_B()): a single whole number of at least the count the
  # interval method and the level require (400 for the percentile interval
  # at 95 percent). Through version 1.0.1 the floor here was 1, and B = 1
  # returned an interval with NA limits. A default B (NULL) is held to the
  # same rule where it is resolved, below.
  if (!is.null(B)) {
    if (ci_supplied && ci_method != "none") {
      .check_boot_B(B, ci_method, conf_level)
    } else if (!is.numeric(B) || length(B) != 1L || !is.finite(B) ||
               B < 1 || B %% 1 != 0) {
      .dmar_refuse_argument(c(
        "'B' must be NULL or a single whole number of at least 1, the ",
        "number of bootstrap replications."), arg = "B")
    }
  }

  # Auxiliary variables carry information only through an imputation
  # model, so supplying them implies multiple imputation, as supplying
  # them to reliability_omega() implies full information maximum
  # likelihood; with a deletion method named explicitly, they are an
  # error rather than silently ignored.
  if (!is.null(aux)) {
    if (missing_supplied && missing != "mi") {
      .dmar_refuse_argument(c(
        "'aux' cannot be combined with missing = \"", missing, "\": ",
        "auxiliary variables carry information only through the ",
        "imputation model. Drop 'aux' or use missing = \"mi\"."), arg = "aux")
    }
    missing <- "mi"
  }
  if (missing != "mi" && (impute_supplied || imputations_supplied)) {
    .dmar_refuse_argument(
      "'impute' and 'imputations' apply only to missing = \"mi\".",
      arg = c("impute", "imputations"))
  }
  if (missing == "mi") {
    if (!is.numeric(imputations) || length(imputations) != 1L ||
        !is.finite(imputations) || imputations < 2 ||
        imputations %% 1 != 0) {
      .dmar_refuse_argument(c(
        "'imputations' must be a single whole number of at least 2: ",
        "the bootstrap's estimate of the fraction of missing ",
        "information compares the imputations within each ",
        "replication, which takes two."), arg = "imputations")
    }
    imputations <- as.integer(imputations)
    if (ci_supplied && ci_method == "bca") {
      .dmar_refuse_argument(c(
        "The BCa interval has not been studied with bootstrap then ",
        "impute, and its jackknife would cost n times 'imputations' ",
        "imputations; use ci_method = \"percentile\" (Schomaker & ",
        "Heumann, 2018)."), arg = "imputations")
    }
  }

  # SPSS user-defined missing values are missing under every treatment
  # (.relia_declared_missing()); through version 1.0.1 a code declared
  # missing was analyzed as a category.
  data <- .relia_declared_missing(data)

  if (missing == "listwise") {
    data <- as.data.frame(data)
    keep <- stats::complete.cases(data)
    # With no complete row there is nothing to fit, and lavaan would
    # print its table of variables and stop with an error that the fit
    # below reports as non-convergence; say what happened and what to
    # use instead.
    if (!any(keep)) .stop_omega_c_no_complete_rows(data)
    all_rows <- data
    data <- data[keep, , drop = FALSE]
    item_data <- data
  } else {
    # Rows with no observed item carry no information about the items
    # and are set aside, unless an auxiliary variable is observed on
    # them (multiple imputation then imputes their items from it), as
    # under full information maximum likelihood in the rest of the
    # family; .relia_resolve_inputs() does that, and validates the
    # auxiliary variables (named columns of `data`, numeric, not
    # constant). The items are ordered categories, checked below, so its
    # checks of continuous items are off (continuous = FALSE).
    inp <- .relia_resolve_inputs(data, S = NULL, N = NULL, allow_S = FALSE,
                                 missing = "fiml", aux = aux,
                                 continuous = FALSE)
    item_data <- inp$data[inp$items]
  }
  # The single-factor model for ordered items is not identified with two
  # items: their one polychoric correlation determines only the product
  # of the two loadings, while categorical omega depends on each, so a
  # two-item fit returns whatever the optimizer's start values lead to.
  # Refuse it, as reliability_omega() refuses two continuous items.
  if (ncol(item_data) < 3L) .stop_omega_c_two_items()
  # Every item as the whole numbers that code its categories: a factor's
  # labels read as numbers, so that listwise and pairwise deletion and
  # multiple imputation rank its categories alike (.omega_c_item_codes()).
  item_data <- .omega_c_item_codes(item_data)
  if (missing != "listwise") inp$data[inp$items] <- item_data

  if (missing != "listwise") {
    unobserved <- colnames(item_data)[colSums(!is.na(item_data)) == 0L]
    if (length(unobserved)) {
      .dmar_refuse_argument(c(
        "No row observes ", paste(unobserved, collapse = " or "),
        ", so neither its thresholds nor its correlations with the ",
        "other items can be estimated. Leave it out of 'data'."), arg = "data")
    }
  }
  # An item with one observed value has no threshold to estimate, and
  # imputation could only repeat that value.
  .check_omega_c_categories(item_data, missing,
                            if (missing == "listwise") all_rows)

  # With every item value observed there is nothing for pairwise
  # deletion or imputation to do: the result is the complete data one,
  # computed as missing = "listwise" computes it, with no imputation.
  if (missing != "listwise" && !anyNA(item_data)) {
    message("No item value is missing",
            if (nrow(item_data) < nrow(as.data.frame(data))) {
              " in the rows with an observed item"
            },
            ", so missing = \"", missing, "\" has nothing to do; the ",
            "result is the complete data one (missing = \"listwise\")",
            # Nothing is imputed, so a replication costs one refit and
            # the listwise default of B applies, not the 1000 of "mi".
            if (missing == "mi" && ci_supplied && ci_method != "none") {
              paste0(", whose bootstrap has B = ",
                     format(if (is.null(B)) 10000 else B,
                            scientific = FALSE),
                     " replications",
                     if (is.null(B)) {
                       paste0(" (its default; each is one refit, with ",
                              "nothing to impute)")
                     })
            },
            ".")
    args <- list(data = item_data, conf_level = conf_level, B = B,
                 seed = seed, denominator = denominator, ncpus = ncpus,
                 progress = progress)
    if (ci_supplied) args$ci_method <- ci_method
    return(do.call(reliability_omega_categorical, args))
  }

  # Each item's categories are scored as consecutive ranks; say so when
  # the codes as given are not one linear rescaling of those ranks shared
  # by every item, since their sum is then another composite.
  .message_omega_c_code_spacing(item_data)

  # B counts bootstrap replications. Each replication of bootstrap then
  # impute costs `imputations` imputations and as many refits, so its
  # default is 1000 rather than the 10000 of the other two treatments.
  if (is.null(B)) {
    B <- if (missing == "mi") 1000 else 10000
    # The default is held to the rule as a supplied B is: 1000 serves the
    # percentile interval through the 95 percent level and not at 99.
    if (ci_method != "none") .check_boot_B(B, ci_method, conf_level)
  }

  N <- nrow(item_data)
  J <- ncol(item_data)
  N_complete <- sum(stats::complete.cases(item_data))

  # Every interval for categorical omega is bootstrap based, and no
  # bootstrap runs unless the user asks for one. The message naming the
  # call that produces the interval waits until the point estimate
  # exists, so a call that is about to error does not first receive
  # advice.
  if (!ci_supplied) ci_method <- "none"

  if (missing == "mi") {
    res <- .omega_c_mi(data = inp$data, items = inp$items, impute = impute,
                       imputations = imputations,
                       denominator = denominator, kind = ci_method,
                       B = B, conf_level = conf_level, seed = seed,
                       ncpus = ncpus, progress = progress)
    # Under the observed denominator a completed data set that leaves a
    # correlation undetermined is set aside, so only the model implied
    # one averages such sets, and only those whose determined
    # correlations identify the loadings (.omega_c_mi_point()).
    if (nrow(res$undetermined) > 0L) {
      one <- nrow(res$undetermined) == 1L
      message("The data do not determine ",
              .undetermined_message(res$undetermined), " in some of the ",
              "completed data sets the estimate averages. The model ",
              "implied denominator does not need ",
              if (one) "it" else "them", ", and in those sets the ",
              "correlations the data do determine identify the ",
              "single-factor loadings; ",
              if (one) "its residual correlation is" else
                "their residual correlations are", " reported as NA.")
    }
    if (!ci_supplied || ci_method == "none") {
      message(
        if (!ci_supplied) {
          paste0("Categorical omega is reported without a confidence ",
                 "interval by default because its interval is bootstrap ",
                 "based. Under missing = \"mi\" request it with ",
                 "ci_method = \"percentile\" (bootstrap then impute); ",
                 "B = 1000 replications, each imputed ", imputations,
                 " times, is the default when you do. ")
        },
        "Without the bootstrap the point estimate averages ",
        if (res$mi$imputations_returned < res$mi$imputations_point) {
          paste0("the ", res$mi$imputations_returned, " of its ",
                 res$mi$imputations_point, " imputations that returned a ",
                 "value")
        } else {
          paste0(res$mi$imputations_point, " imputations")
        },
        "; the number von Hippel's (2020) rule asks for needs the ",
        "fraction of missing information, which the bootstrap estimates.")
    }
    if (length(res$notes)) {
      # Classed as the bootstrap engine classes its warning: an interval
      # withheld by the share or the count rule, or replications set
      # aside with the interval formed; anything else is a
      # dmar_bootstrap_note.
      # It names the function, as every warning of the engine does.
      text <- paste0("reliability_omega_categorical(): ",
                     paste(res$notes, collapse = ". "), ".")
      bs <- res$bootstrap
      fields <- if (!is.null(bs$bootstrap_status)) list(
        fun = "reliability_omega_categorical",
        B = as.integer(bs$attempted), B_used = bs$B_used,
        n_dropped = as.integer(bs$attempted) - bs$B_used,
        B_required = bs$required,
        share_dropped = (bs$attempted - bs$B_used) / bs$attempted,
        caution = isTRUE(bs$caution), status = bs$bootstrap_status)
      if (!is.null(fields) && startsWith(bs$bootstrap_status, "insufficient")) {
        do.call(.dmar_warn, c(list(text, "bootstrap_insufficient"), fields))
      } else if (!is.null(fields) && isTRUE(bs$failed > 0L)) {
        do.call(.dmar_warn,
                c(list(text, "bootstrap_replications_dropped"), fields))
      } else {
        do.call(.dmar_warn, c(list(text, "bootstrap_note"), fields))
      }
    }
    if (identical(res$bootstrap$bootstrap_status, "exploratory")) {
      bs <- res$bootstrap
      marked <- .dmar_boot_settle(
        bs$attempted, bs$B_used, kind = .boot_method_kind(ci_method),
        tail = (1 - conf_level) / 2, conf_level = conf_level,
        fun = "reliability_omega_categorical",
        quantity = "categorical omega",
        dropped_because = "returned no value", signal = FALSE)
      .dmar_warn(.dmar_boot_message(marked, "exploratory"),
                 "bootstrap_exploratory", B = marked$B,
                 B_used = marked$B_used, n_dropped = marked$n_dropped,
                 B_required = marked$B_required,
                 share_dropped = marked$share_dropped,
                 caution = marked$caution, status = marked$status)
    }
    out <- .relia_result(
      estimate = res$estimate, se = res$se,
      lower = res$lower, upper = res$upper,
      conf_level = conf_level, N = N, J = J,
      coefficient = "omega_categorical",
      ci_method = res$ci_method,
      bootstrap = res$bootstrap,
      B = if (ci_method != "none") B else NA_integer_,
      N_complete = N_complete, missing = "mi", aux = inp$aux
    )
    attr(out, "denominator") <- denominator
    attr(out, "imputation") <- res$mi
    if (!is.null(res$residual_cor)) {
      attr(out, "residual_cor") <- res$residual_cor
    }
    return(out)
  }

  if (missing == "pairwise") {
    # A pair of items never observed together has no polychoric
    # correlation under pairwise deletion.
    joint <- crossprod(!is.na(as.matrix(item_data)))
    never <- which(joint == 0 & upper.tri(joint), arr.ind = TRUE)
    if (nrow(never) > 0L) {
      nm <- colnames(item_data)
      .dmar_refuse_argument(c(
        "No row observes both ", nm[never[1L, 1L]], " and ",
        nm[never[1L, 2L]], ", so pairwise deletion has no data for ",
        "their polychoric correlation. Use missing = \"mi\", whose ",
        "imputation model relates them through the other items."), arg = NULL)
    }
  }

  # One fit supplies the point estimate and the residual polychoric
  # correlations; the bootstrap refits through .omega_c_from_data() with
  # the same denominator and missing data treatment on every
  # replication.
  # A constant item has no estimable thresholds. Refuse it by name
  # before lavaan prints its variable table, as for continuous resamples.
  point_fn <- .relia_checked_point(function(d) {
    .omega_c_from_data(d, denominator = denominator, missing = missing)
  }, missing = missing)
  parts <- .omega_c_parts(item_data, missing = missing)
  # No fit: lavaan's own message when it stopped with an error, which
  # through version 1.0.1 was reported as a model that did not converge,
  # and otherwise that the model did not converge.
  if (!isTRUE(parts$converged)) {
    .dmar_refuse(
      .relia_fit_failure(parts, model = "single-factor categorical CFA",
                         missing = missing), "computation_failed")
  }
  undetermined <- parts$undetermined
  # Correlations the data do not determine carry no weight in the fit,
  # and when the determined ones do not connect every item to the others
  # and include a cycle of an odd number of items the loadings are not
  # identified (.omega_c_identification()): categorical omega under
  # either denominator would be whatever the start values led to, which
  # through the first versions of 1.1.0 the model implied denominator
  # returned. This comes first, since the loadings an improper verdict
  # reads are as arbitrary as the rest.
  if (isFALSE(parts$identified)) {
    .dmar_refuse(
      .omega_c_unidentified_message(parts$identification, undetermined),
      "computation_failed")
  }
  # An improper solution is reported before anything about undetermined
  # pairs, so no advice to switch denominators precedes an error that
  # either denominator would meet.
  if (length(parts$improper)) {
    .dmar_refuse(c(
      .improper_message(parts),
      if (nrow(undetermined) > 0L) {
        paste0(" The data also do not determine ",
               .undetermined_message(undetermined), ".")
      }), "computation_failed")
  }
  if (nrow(undetermined) > 0L) {
    what <- .undetermined_message(undetermined)
    one <- nrow(undetermined) == 1L
    table_has <- if (one) "the bivariate table has" else
      "the bivariate tables have"
    if (denominator == "observed") {
      # The model implied denominator is offered only here, where the
      # determined correlations identify the loadings (checked above).
      .dmar_refuse_argument(c(
        "The data do not determine ", what, ": ", table_has, " empty ",
        "cells (two items skewed in opposite directions produce this), ",
        "so the standard ", if (one) "error is" else "errors are",
        " unbounded. The observed denominator needs every unrestricted ",
        "polychoric correlation. Use denominator = \"model_implied\" ",
        "(Green & Yang, 2009), which takes the ",
        if (one) "correlation" else "correlations",
        " the single-factor model implies for ",
        if (one) "that pair" else "those pairs",
        "; the correlations the data do determine identify the loadings ",
        "here."), arg = NULL)
    }
    message("The data do not determine ", what, " (", table_has, " empty ",
            "cells). The model implied denominator does not need ",
            if (one) "it" else "them", ", and the correlations the data do ",
            "determine identify the single-factor loadings; ",
            if (one) "its residual correlation is" else
              "their residual correlations are", " reported as NA.")
  }
  rho_den <- if (denominator == "model_implied") {
    parts$rho_model
  } else {
    parts$rho_sat
  }
  estimate <- .gy_omega_from_components(parts$thresh, parts$rho_num,
                                        rho_den)
  if (!is.finite(estimate)) {
    .dmar_refuse(
      "The single-factor categorical CFA model did not converge.",
      "not_converged")
  }
  resid <- parts$rho_sat - parts$rho_model
  for (k in seq_len(nrow(undetermined))) {
    resid[undetermined[k, 1], undetermined[k, 2]] <- NA_real_
    resid[undetermined[k, 2], undetermined[k, 1]] <- NA_real_
  }
  residual_cor <- .residual_cor_summary(resid, df = parts$df)
  if (!ci_supplied) {
    message("Categorical omega is reported without a confidence interval ",
            "by default because its interval is bootstrap based. Request ",
            "it with ci_method = \"bca\" (the recommended method) or ",
            "\"percentile\"; B = 10000 replications is the default when ",
            "you do.")
  }

  ci <- switch(
    ci_method,
    none = list(se = NA_real_, lower = NA_real_, upper = NA_real_),
    bootstrap_se = ,
    bootstrap_se_logistic = ,
    percentile = ,
    # Categorical omega is not defined for an improper solution, so an
    # improper replication has no estimate to keep and is set aside
    # under either rule; "drop" records that in the bootstrap attribute.
    # No interval is formed from fewer than half of the replications
    # (see .bootstrap_ci()), the criterion the point
    # estimate under multiple imputation applies to its completed data
    # sets.
    bca = .bootstrap_ci(
      fun = "reliability_omega_categorical",
      data = item_data,
      point_fn = point_fn,
      B = B, conf_level = conf_level, kind = ci_method, seed = seed,
      ncpus = ncpus, improper = "drop", progress = progress,
      undetermined_note = .omega_c_undetermined_note(denominator)
    )
  )

  out <- .relia_result(
    estimate = estimate, se = ci$se,
    lower = ci$lower, upper = ci$upper,
    conf_level = conf_level, N = N, J = J,
    coefficient = "omega_categorical",
    ci_method = ci$ci_method %||% ci_method,
    bootstrap = ci$bootstrap,
    B = if (ci_method != "none") B else NA_integer_,
    N_complete = N_complete, missing = missing
  )
  attr(out, "denominator") <- denominator
  if (!is.null(residual_cor)) attr(out, "residual_cor") <- residual_cor
  out
}


# .stop_omega_c_no_complete_rows(data)
#
# The refusal when listwise deletion leaves fewer than two rows: it
# names an item no row observes (for which no treatment helps), and
# otherwise the treatments that use the incomplete rows, offering
# pairwise deletion only when every pair of items is observed together
# somewhere, since without that it has no data for some polychoric
# correlation either.
.stop_omega_c_no_complete_rows <- function(data) {
  data <- as.data.frame(data)
  obs <- !is.na(as.matrix(data))
  nm <- colnames(data)
  unobserved <- nm[colSums(obs) == 0L]
  if (length(unobserved)) {
    .dmar_refuse_argument(c(
      "No row observes ", paste(unobserved, collapse = " or "),
      ", so neither its thresholds nor its correlations with the ",
      "other items can be estimated. Leave it out of 'data'."), arg = "data")
  }
  joint <- crossprod(obs)
  never <- which(joint == 0 & upper.tri(joint), arr.ind = TRUE)
  remedy <- if (nrow(never) > 0L) {
    paste0("Use missing = \"mi\": no row observes both ",
           nm[never[1L, 1L]], " and ", nm[never[1L, 2L]], ", so ",
           "pairwise deletion has no data for their polychoric ",
           "correlation either.")
  } else {
    paste0("Use missing = \"pairwise\", which assumes values missing ",
           "completely at random as listwise deletion does, or missing = ",
           "\"mi\", which assumes them missing at random.")
  }
  n_complete <- sum(stats::complete.cases(data))
  .dmar_refuse_argument(c(
    if (n_complete == 0L) "No row" else "Only one row", " of 'data' ",
    "is observed on every item, so listwise deletion (the default) ",
    "leaves nothing to analyze. ", remedy), arg = "data")
}


# .check_omega_c_categories(item_data, missing, all_rows = NULL)
#
# Stops when an item takes a single value in the rows analyzed (the
# complete rows under listwise deletion, the observed values otherwise):
# it has no threshold to estimate and no variance to share with the
# factor, lavaan stops on it with an error the caller would report as
# non-convergence, and imputation could only repeat the one value. With
# fewer than two complete rows every item does, which
# .stop_omega_c_no_complete_rows() says instead. Under listwise deletion
# the message notes when the item varies in the rows deletion removed.
#
# Inputs:
#   item_data : the item columns of the rows analyzed.
#   missing   : the treatment.
#   all_rows  : under listwise deletion, the data before deletion.
.check_omega_c_categories <- function(item_data, missing, all_rows = NULL) {
  if (missing == "listwise" && nrow(item_data) < 2L) {
    .stop_omega_c_no_complete_rows(all_rows)
  }
  n_values <- vapply(item_data, function(x) length(unique(x[!is.na(x)])),
                     integer(1L))
  single <- names(item_data)[n_values < 2L]
  if (!length(single)) return(invisible(TRUE))
  who <- if (length(single) == 1L) {
    paste(single, "takes")
  } else {
    paste(paste(single, collapse = " and "), "each take")
  }
  if (missing == "listwise") {
    varies <- !is.null(all_rows) && any(vapply(single, function(j) {
      length(unique(all_rows[[j]][!is.na(all_rows[[j]])])) > 1L
    }, logical(1L)))
    .dmar_refuse_argument(c(
      who, " a single value in the ", nrow(item_data), " complete ",
      "rows, so ", if (length(single) == 1L) "its" else "their",
      " thresholds cannot be estimated.",
      if (varies) {
        paste0(" Listwise deletion (the default) analyzes only the ",
               "complete rows; missing = \"pairwise\" or \"mi\" uses ",
               "the incomplete rows as well.")
      } else {
        " Leave it out of 'data'."
      }), arg = NULL)
  }
  .dmar_refuse_argument(c(
    who, " a single observed value, so ",
    if (length(single) == 1L) "its" else "their", " thresholds cannot ",
    "be estimated",
    if (missing == "mi") ", and imputation could only repeat it",
    ". Leave it out of 'data'."), arg = "data")
}


# .omega_c_item_codes(item_data)
#
# The items as the whole numbers that code their categories. A numeric
# column, a labeled one (haven) included, is read through its plain
# values and must hold whole numbers. A factor, or a character column,
# is read through its labels, which must name whole numbers: a Likert
# factor labeled 0 to 4, or one whose lowest category no respondent
# chose (levels "2" to "5"), is the item its numbers code, and is
# returned as those numbers, so .omega_c_parts() orders its categories
# by them under listwise and pairwise deletion, as multiple imputation
# does (.mi_plain_columns()). Through version 1.0.1 a factor passed only
# when each label equaled its level's index (1, 2, ..., K), and any
# other factor, its labels whole numbers or not, was refused with
# "column j contains non-integer values". An ordered factor whose levels
# are not in the order of the numbers they name is refused: its order
# and its codes disagree, and the function would have to choose one. A
# factor whose labels are words (strongly disagree, ...) is refused with
# a message that says how to give its codes.
#
# Inputs:
#   item_data : data frame of the item columns.
# Returns:
#   item_data, with every factor or character column replaced by the
#   integers its labels name (NA kept).
.omega_c_item_codes <- function(item_data) {
  shown <- function(v) {
    v <- utils::head(v, 5L)
    paste0("\"", v, "\"", collapse = ", ")
  }
  for (j in seq_len(ncol(item_data))) {
    x <- item_data[[j]]
    nm <- names(item_data)[j]
    where <- if (is.null(nm) || !nzchar(nm)) {
      paste("column", j)
    } else {
      paste0(nm, " (column ", j, ")")
    }
    if (is.factor(x) || is.character(x)) {
      lab <- if (is.factor(x)) levels(x) else sort(unique(x[!is.na(x)]))
      num <- suppressWarnings(as.numeric(lab))
      if (anyNA(num) || any(num != round(num))) {
        .dmar_refuse_argument(c(
          "Categorical omega reads each item's categories from whole ",
          "number codes, and ", where, " is a ",
          if (is.factor(x)) "factor" else "character column",
          " whose labels are not all whole numbers (", shown(lab), "). ",
          "Give its categories as integer codes in their order, such as ",
          "as.integer(x) for a factor whose levels are in order."), arg = NULL)
      }
      if (is.ordered(x) && is.unsorted(num, strictly = TRUE)) {
        .dmar_refuse_argument(c(
          where, " is an ordered factor whose levels (", shown(lab),
          ") are not in the order of the numbers they name; categorical ",
          "omega orders an item's categories by their codes. Give the ",
          "codes in the order of the categories, or relevel the factor."),
          arg = NULL)
      }
      item_data[[j]] <- as.integer(num[match(as.character(x), lab)])
      next
    }
    v <- as.vector(unclass(x))
    if (!is.numeric(v) && !is.logical(v) ||
        !all(v == round(v), na.rm = TRUE)) {
      .dmar_refuse_argument(c(
        "Categorical omega requires items coded as integer category ",
        "values; ", where, " contains non-integer values."), arg = NULL)
    }
  }
  item_data
}


# .message_omega_c_code_spacing(item_data)
#
# Categorical omega is the reliability of the sum of category ranks:
# .omega_c_parts() makes each item an ordered factor of its observed
# codes, so its categories are scored 0, 1, ..., C - 1 in the order of
# the codes (Kelley & Pornprasertmanit, 2016, Eq. 16, X = c with c = 0,
# 1, ..., C - 1; Green & Yang, 2009, score the categories the same way),
# and MBESS and semTools compute the same coefficient. The sum of the
# codes as given is a linear function of that sum, and so has the same
# reliability, exactly when every item's observed codes are equally
# spaced with one step shared by all items (0, 2, 4 on every item, say).
# Otherwise (codes such as -2, -1, 1, 2; items stepping by 2 beside
# items stepping by 1; a category between two used ones that no
# respondent chose) the coded sum is another composite, and this message
# says so; the coefficient is not changed. A resample of the bootstrap
# that loses an interior category is not checked: its effect on the
# coefficient is of the order of 0.001.
#
# Inputs:
#   item_data : the item columns of the rows analyzed (integer codes).
# Returns:
#   invisible(NULL); gives a message when the codes are not one shared
#   linear rescaling of the ranks.
.message_omega_c_code_spacing <- function(item_data) {
  codes <- lapply(item_data, function(x) {
    x <- if (is.factor(x)) as.character(x) else as.vector(unclass(x))
    sort(unique(as.numeric(x[!is.na(x)])))
  })
  steps <- lapply(codes, diff)
  uneven <- names(codes)[vapply(steps, function(s) {
    length(unique(s)) > 1L
  }, logical(1L))]
  step <- vapply(steps, function(s) {
    if (length(unique(s)) == 1L) s[1L] else NA_real_
  }, numeric(1L))
  even_steps <- step[!is.na(step)]
  differ <- length(unique(even_steps)) > 1L
  if (!length(uneven) && !differ) return(invisible(NULL))
  shown <- function(j) {
    v <- codes[[j]]
    txt <- format(v, trim = TRUE, scientific = FALSE)
    if (length(txt) > 8L) txt <- c(txt[1:7], "...")
    paste0(j, " is coded ", paste(txt, collapse = ", "))
  }
  detail <- character()
  if (length(uneven)) {
    detail <- c(detail, paste(vapply(utils::head(uneven, 3L), shown,
                                     character(1L)), collapse = "; "))
  }
  if (differ) {
    by_step <- split(names(even_steps),
                     factor(even_steps, levels = unique(even_steps)))
    detail <- c(detail, paste(vapply(names(by_step), function(k) {
      paste0(paste(by_step[[k]], collapse = " and "), " step by ", k)
    }, character(1L)), collapse = ", "))
  }
  message("Categorical omega scores the categories of each item as ",
          "consecutive ranks, 0, 1, ..., C - 1 in the order of its codes ",
          "(Green & Yang, 2009; Kelley & Pornprasertmanit, 2016), so it is ",
          "the reliability of the sum of those ranks. These codes are not ",
          "one rescaling of the ranks shared by every item (",
          paste(detail, collapse = "; "), "), so the sum of the codes as ",
          "given is a different composite, whose reliability can differ.")
  invisible(NULL)
}
