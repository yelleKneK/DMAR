#' Controlled Test-Market Experiment (Bryant & Bruvold, 1980)
#'
#' The controlled test-market experiment of Bryant and Bruvold (1980),
#' used to illustrate multiple-comparison procedures in the analysis of
#' covariance (ANCOVA) when the covariate is \emph{random}. A company
#' compared \eqn{k = 6} marketing strategies (\dQuote{panels}) for a brand,
#' randomly assigning them to retail outlets within \eqn{s = 4} blocks of
#' outlets that were homogeneous in size, locality, and ownership (a
#' randomized complete block design, one outlet per panel-by-block cell).
#' During the experiment a concomitant variable -- the remaining category
#' movement in each outlet -- becomes available; it cannot be controlled by
#' the experimenter and is best modeled as a random covariate. Adjusting
#' brand movement for this covariate sharply reduces unexplained error,
#' permitting far finer comparison of the panels than the raw outcome allows.
#'
#' @format A data frame with 24 observations (6 panels \eqn{\times} 4 blocks)
#' on 4 variables.
#' \describe{
#'   \item{\code{panel}}{Factor with levels \code{1}--\code{6}: the marketing
#'     strategy (treatment) randomly assigned to the outlet. Different panels
#'     entail different methods of packaging, displaying, or pricing.}
#'   \item{\code{block}}{Factor with levels \code{1}--\code{4}: the block of
#'     retail outlets, grouped to be homogeneous in size, locality, ownership,
#'     and other considerations that influence brand movement.}
#'   \item{\code{brand_movement}}{Test-brand movement during the test period,
#'     in hundreds of statistical cases. The dependent variable (\eqn{y} in
#'     the source).}
#'   \item{\code{category_movement}}{Remaining category movement, the random
#'     concomitant variable (covariate; \eqn{x} in the source). It is not
#'     identically distributed across blocks, which is precisely the setting
#'     Bryant and Bruvold's grouped-covariate extension was designed for.}
#' }
#'
#' @details
#' \strong{Why it is a benchmark for ANCOVA multiple comparisons.} The
#' model fitted by Bryant and Bruvold (their Eq. 3.1) is a randomized-block
#' ANCOVA,
#' \deqn{y_{ij} = \theta_i + \beta_j + (x_{ij} - \delta_j)\, u + e_{ij},}
#' with \eqn{\theta_i} the \eqn{i}th panel (adjusted) mean, \eqn{\beta_j} the
#' \eqn{j}th block effect, and \eqn{u} the within-cell covariate slope. The
#' point of the example is that the studentized range of the adjusted panel
#' means does \emph{not} follow the ordinary Tukey distribution, because the
#' covariate is random and its adjustment must be estimated; the correct
#' reference distribution is the Bryant--Paulson generalized studentized
#' range (\code{\link{bryant_paulson}}).
#'
#' \strong{Reproducible quantities.} Fitting
#' \code{lm(brand_movement ~ panel + block + category_movement)} gives a
#' covariate slope of \eqn{0.4079} and an error mean square of
#' \eqn{0.01326} on \eqn{\nu = 14} degrees of freedom, with adjusted panel
#' means \eqn{3.595, 3.619, 4.102, 4.515, 4.618, 4.876} -- exactly the values
#' reported in the paper. With \eqn{q_{.05;\,1,6,14} = 4.83}
#' (\code{\link{qbryant_paulson}}), every pairwise simultaneous 95\%
#' interval is a difference of adjusted panel means plus or minus
#' \eqn{0.278}, so two panels differ at the simultaneous 95\% level
#' exactly when their adjusted means are more than \eqn{0.278} apart.
#' Had the covariate not been measured, the error mean square
#' would have been \eqn{0.2368} -- roughly eighteen times larger -- and the
#' intervals about four times wider. See \code{data-raw/test_market.R} for
#' the construction script and its verification checks.
#'
#' @author Ken Kelley
#'
#' @source
#' Bryant, J. L., & Bruvold, N. T. (1980). Multiple comparison procedures in
#' the analysis of covariance. \emph{Journal of the American Statistical
#' Association, 75}(372), 874--880 (Table 1). \doi{10.2307/2287175}
#'
#' @references
#' Bryant, J. L., & Paulson, A. S. (1976). An extension of Tukey's method of
#'   multiple comparisons to experimental designs with random concomitant
#'   variables. \emph{Biometrika, 63}, 631--638.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 9.)
#'
#' @seealso \code{\link{ci_c_ancova_bp}} for the simultaneous intervals this
#'   data set illustrates, \code{\link{bryant_paulson}} for the critical
#'   values, and \code{\link{ancova}} for an ANCOVA fit.
#'
#' @examples
#' data(test_market)
#' str(test_market)
#'
#' # Reproduce the published ANCOVA (slope 0.4079, error MS 0.01326, df 14).
#' fit <- lm(brand_movement ~ panel + block + category_movement,
#'           data = test_market)
#' coef(fit)["category_movement"]
#' sum(residuals(fit)^2) / fit$df.residual
#'
#' # Adjusted panel means at the covariate grand mean.
#' xbar <- mean(test_market$category_movement)
#' adj <- vapply(levels(test_market$panel), function(p) {
#'   nd <- data.frame(panel = factor(p, levels = levels(test_market$panel)),
#'                    block = factor(1:4, levels = levels(test_market$block)),
#'                    category_movement = xbar)
#'   mean(predict(fit, nd))
#' }, numeric(1))
#' adj  # 3.595 3.619 4.102 4.515 4.618 4.876
#'
#' # Bryant-Paulson simultaneous 95% intervals (s = 4 blocks => n = 4,
#' # df 14), every one of them the difference plus or minus 0.278. Not
#' # run here, because the Bryant-Paulson critical value is obtained by
#' # inverting an integral with uniroot, which takes about half a
#' # second; the call is:
#' # ci_c_ancova_bp(adj_means = adj, s_ancova = sqrt(0.01326),
#' #                n = 4, num_covariates = 1, df = 14)
#' @keywords datasets
"test_market"
