#' Maximum Likelihood Covariance Matrix and Means of Incomplete Data
#'
#' @description
#' For continuous variables some of whose values are missing,
#' \code{cov_fiml()} returns the maximum likelihood estimates of their
#' covariance matrix and mean vector under the multivariate normal model,
#' computed from every observed value rather than from the complete cases
#' alone. These are the estimates of the saturated (unrestricted) model
#' fitted by full information maximum likelihood (FIML), the model in which
#' every variance, covariance, and mean is free. They are computed by the
#' expectation maximization (EM) algorithm of Dempster, Laird, and Rubin
#' (1977) and are often called the EM covariance matrix and the EM means
#' (Enders & Peugh, 2004). The function requires \pkg{lavaan}, whose EM
#' implementation does the computing.
#'
#' The usual alternatives each give something up. Listwise deletion keeps
#' only the complete cases, and its estimates are consistent in general
#' only when values are missing completely at random (MCAR). Pairwise
#' deletion estimates each covariance from a different subset of cases, so
#' the matrix need not be positive definite and has no single sample size
#' (Marsh, 1998). The EM estimates use every case and are consistent when
#' values are missing at random (MAR) given the variables in the model
#' (Schafer & Graham, 2002). They are therefore the matrix to hand to a
#' method that starts from a covariance or correlation matrix, such as an
#' exploratory factor analysis or a parallel analysis of incomplete items.
#'
#' @param data A \code{data.frame} or numeric matrix of continuous
#'   variables, with missing values coded \code{NA}. Every column that
#'   \code{aux} does not name is an analysis variable, a row and column of
#'   the returned \code{S}. Every column must be numeric; a factor,
#'   character, or logical column is refused, as is a column with fewer
#'   than two distinct observed values.
#' @param aux Optional character vector naming auxiliary columns of
#'   \code{data}: variables that are not wanted in \code{S} but carry
#'   information about the missing values, because they predict which
#'   values are missing, the missing values themselves, or both. They enter
#'   the saturated model with the analysis variables, so their information
#'   reaches \code{S} and \code{mu}, and are then left out of what is
#'   returned (the saturated correlates treatment of Graham, 2003, the one
#'   the reliability family's \code{aux} argument uses).
#'
#' @details
#' \strong{What is estimated.} With complete data the estimates are the
#' sample moments with the maximum likelihood divisor: \code{S} equals
#' \code{cov(data) * (N - 1) / N} and \code{mu} equals
#' \code{colMeans(data)}. With missing values, the EM algorithm alternates
#' two steps until the estimates stop changing. The expectation step fills
#' in the sums and cross products each case would contribute, replacing a
#' missing value by its regression prediction from the observed values of
#' its case under the current estimates and adding the residual
#' covariance to the cross products; the maximization step recomputes the
#' moments from those sums. No iteration lowers the likelihood of the
#' observed data. Meeting the stopping criterion establishes numerical
#' stability of the iterations, not a finite global maximum or uniqueness.
#' The estimates are those of the unrestricted (h1) model of
#' every FIML fit in \pkg{lavaan}:
#' \code{lavaan::lavInspect(fit, "sampstat.h1")} returns them for a model
#' fitted with \code{missing = "ml"} to the same variables.
#' \code{cov_fiml()} computes them with \pkg{lavaan}'s EM
#' (\code{lavaan::lavCor()} with \code{missing = "ml"}) and its default
#' convergence criterion, which stops the iterations when no mean or
#' covariance changes by more than \eqn{10^{-5}}. That is a criterion on
#' the absolute change, and in the units of variables whose standard
#' deviations are near .001 it is met after an iteration or two, far from
#' the maximum. \code{cov_fiml()} therefore runs the algorithm on the
#' variables centered at their means and divided by their standard
#' deviations, each computed from the variable's observed values, and
#' transforms the estimates back, so that they do not depend on the units:
#' multiplying a variable by a constant multiplies its mean, and its row
#' and column of \code{S}, by that constant. They agree with
#' \pkg{lavaan}'s own estimates of the same moments to the precision of
#' the criterion, about seven significant digits. When the algorithm
#' stops at its iteration limit instead, \code{converged} is
#' \code{FALSE}, a warning says so, and \code{S} and \code{mu} are the
#' last iterate rather than the maximum likelihood estimates.
#'
#' \strong{When the estimates are consistent.} The estimates are
#' consistent when the values are MAR given the variables in the saturated
#' model: the probability that a value is missing may depend on the
#' observed values of the variables included, but not on the missing
#' values themselves (Schafer & Graham, 2002). MAR is defined relative to
#' the variables analyzed (Enders, 2001). If whether an item is missing
#' depends on a variable that is left out, MAR need not hold for the items
#' alone, and naming that variable in \code{aux} brings it back into the
#' model without making it part of \code{S}. Including such variables
#' (the inclusive strategy) costs little and can remove a bias that
#' leaving them out would leave in place (Collins, Schafer, & Kam, 2001).
#' Under missingness that depends on the missing values themselves (not
#' at random), no choice of auxiliaries guarantees consistency. Normality
#' enters as the working likelihood: Yuan, Yang-Wallentin, and Bentler
#' (2012) summarize results showing that the normal theory estimates can
#' remain consistent when the population distribution is unknown, and
#' their simulations found the normal theory variance estimates biased at
#' smaller sample sizes when the data were not normal, a bias that shrank
#' as the sample grew.
#'
#' \strong{Cases used.} A case on which every variable, analysis and
#' auxiliary, is missing carries no information and is set aside, as
#' \pkg{lavaan} would set it aside. Every other case is used, including a
#' case whose analysis variables are all missing but whose auxiliaries are
#' observed. Under MAR given an auxiliary, such cases carry the auxiliary
#' values of exactly the people whose analysis variables are missing,
#' which is how the model learns how those people differ, and setting
#' them aside would pull the estimates back toward those of the cases that
#' responded. \code{N} counts the cases used and \code{N_complete} the
#' cases complete on the analysis variables, all that listwise deletion
#' would keep.
#'
#' \strong{Refusals.} Two variables that are never observed together have
#' a covariance the data say nothing about: the likelihood is flat in
#' that direction, and an EM estimate of it reflects the starting values
#' and the other covariances rather than the data. The same holds, in
#' substance, for two variables observed together only in rows where one
#' of them takes a single value (or in one row only): those rows say
#' nothing about how the two vary together, and the estimate of their
#' covariance rests on the form of the normal model, or is not determined
#' at all when one of them is observed nowhere else. The function refuses
#' either kind of pair and names it; \code{\link{missing_coverage}}
#' reports the number of cases behind every pair.
#'
#' \strong{Using the matrix downstream.} \code{N} is the number of cases
#' used, not an effective sample size, and the matrix carries less
#' information than \code{N} complete cases would. A method that takes a
#' covariance matrix and a sample size and treats them as complete data
#' gets consistent point estimates from \code{S}, but its standard errors
#' and test statistics are too optimistic, because they ignore the
#' information the missing values took away (Savalei & Bentler, 2009).
#' Enders and Peugh (2004) studied which sample size to give such a
#' method when fitting a confirmatory factor model: the smallest number of
#' cases behind any covariance gave Type I error rates near the nominal
#' level when every indicator had missing values (the smallest number
#' behind any variance did when only some indicators had them), and the
#' harmonic mean of the numbers behind the variances gave confidence
#' intervals with coverage near the nominal level.
#' \code{\link{missing_coverage}} reports the number of cases behind every
#' variance and covariance. The two-stage approach of Savalei and Bentler
#' (2009) corrects the standard errors and the test statistic of a model
#' fitted to the EM moments (\pkg{lavaan}'s \code{missing = "two.stage"}),
#' and a model fitted to the raw data by FIML needs no correction. The
#' fraction of missing information
#' (\code{\link{fraction_missing_information}}) says what the missing
#' values cost each parameter of such a fit. For point estimates, such as
#' the eigenvalues of a factor analysis or of a parallel analysis, the EM
#' matrix is the consistent basis.
#'
#' \strong{Standard errors.} The function returns point estimates only.
#' The sampling variability that matters is that of whatever is computed
#' from the matrix, and the analysis that computes it is where its interval
#' belongs.
#'
#' @return A list with the same six elements whatever the data. A matrix
#'   is what every use of the result takes (\code{stats::cov2cor()},
#'   \code{stats::factanal(covmat = )}, a function's \code{S} argument), so
#'   the result holds the matrix itself, as \code{\link{cov_sem}} and
#'   \code{\link{covmat_from_cfa}} return theirs, rather than a table of
#'   its elements:
#'   \describe{
#'     \item{\code{S}}{The maximum likelihood covariance matrix of the
#'       analysis variables (divisor \code{N}), with their names as row
#'       and column names, in the order of the columns of \code{data}.}
#'     \item{\code{mu}}{The maximum likelihood mean vector of the analysis
#'       variables, named, in the same order.}
#'     \item{\code{N}}{The number of cases used: those with at least one
#'       observed value among the analysis variables and the
#'       auxiliaries.}
#'     \item{\code{N_complete}}{The number of cases observed on every
#'       analysis variable.}
#'     \item{\code{patterns}}{The number of distinct patterns of observed
#'       and missing values among the cases used, over every variable in
#'       the saturated model (the analysis variables and the auxiliaries);
#'       \code{\link{missing_patterns}} tabulates them.}
#'     \item{\code{converged}}{\code{TRUE} when the EM algorithm met its
#'       convergence criterion.}
#'   }
#'
#' @references
#' Collins, L. M., Schafer, J. L., & Kam, C.-M. (2001). A comparison of
#'   inclusive and restrictive strategies in modern missing data
#'   procedures. \emph{Psychological Methods, 6}(4), 330--351.
#'   \doi{10.1037/1082-989X.6.4.330}
#'
#' Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). Maximum likelihood
#'   from incomplete data via the EM algorithm. \emph{Journal of the Royal
#'   Statistical Society: Series B (Methodological), 39}(1), 1--38.
#'
#' Enders, C. K. (2001). A primer on maximum likelihood algorithms
#'   available for use with missing data. \emph{Structural Equation
#'   Modeling, 8}(1), 128--141.
#'
#' Enders, C. K., & Peugh, J. L. (2004). Using an EM covariance matrix to
#'   estimate structural equation models with missing data: Choosing an
#'   adjusted sample size to improve the accuracy of inferences.
#'   \emph{Structural Equation Modeling, 11}(1), 1--19.
#'
#' Graham, J. W. (2003). Adding missing-data-relevant variables to
#'   FIML-based structural equation models. \emph{Structural Equation
#'   Modeling, 10}(1), 80--100. \doi{10.1207/S15328007SEM1001_4}
#'
#' Marsh, H. W. (1998). Pairwise deletion for missing data in structural
#'   equation models: Nonpositive definite matrices, parameter estimates,
#'   goodness of fit, and adjusted sample sizes. \emph{Structural Equation
#'   Modeling, 5}(1), 22--36. \doi{10.1080/10705519809540087}
#'
#' Savalei, V., & Bentler, P. M. (2009). A two-stage approach to missing
#'   data: Theory and application to auxiliary variables. \emph{Structural
#'   Equation Modeling, 16}(3), 477--497. \doi{10.1080/10705510903008238}
#'
#' Schafer, J. L., & Graham, J. W. (2002). Missing data: Our view of the
#'   state of the art. \emph{Psychological Methods, 7}(2), 147--177.
#'   \doi{10.1037/1082-989X.7.2.147}
#'
#' Yuan, K.-H., Yang-Wallentin, F., & Bentler, P. M. (2012). ML versus MI
#'   for missing data with violation of distribution conditions.
#'   \emph{Sociological Methods & Research, 41}(4), 598--629.
#'   \doi{10.1177/0049124112460373}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{missing_patterns}} and
#'   \code{\link{missing_coverage}} to describe the missing values first;
#'   \code{\link{fraction_missing_information}} for what they cost the
#'   parameters of a FIML fit; \code{\link{reliability_alpha}} and
#'   \code{\link{reliability_omega}}, whose \code{missing = "fiml"} paths
#'   also bring auxiliary variables into a saturated model.
#'
#' @family missing data
#'
#' @keywords multivariate
#'
#' @examples
#' # The six spatial tests of the Holzinger and Swineford data. The revised
#' # paper form board (t25) and flags (t26) tests were given only at
#' # Grant-White, so the 156 Pasteur pupils are missing both of them.
#' spatial <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
#'              "t4_lozenges", "t25_paper_form_board_r", "t26_flags")
#' em <- cov_fiml(holzinger_swineford[spatial])
#' em$N
#' em$N_complete
#' em$patterns
#' round(stats::cov2cor(em$S), 2)
#'
#' # Listwise deletion would keep the 145 Grant-White pupils alone and
#' # estimate every correlation from them.
#' round(stats::cor(holzinger_swineford[spatial], use = "complete.obs"), 2)
#'
#' # Values missing at random given a variable that is not analyzed: word
#' # meaning (t9) goes missing more often for pupils who scored high on
#' # general information (t5). The two visual tests say little about those
#' # pupils, so without t5 the mean of t9 stays near that of the pupils who
#' # took it; named as an auxiliary, t5 recovers most of the difference.
#' set.seed(113)
#' d <- holzinger_swineford[c("t1_visual_perception", "t4_lozenges",
#'                            "t9_word_meaning", "t5_general_information")]
#' z <- as.numeric(scale(d$t5_general_information))
#' p_missing <- stats::plogis(-0.5 + 2 * z)
#' d$t9_word_meaning[stats::runif(nrow(d)) < p_missing] <- NA
#' without_aux <- cov_fiml(d[1:3])
#' with_aux <- cov_fiml(d, aux = "t5_general_information")
#' rbind(without_aux = without_aux$mu, with_aux = with_aux$mu,
#'       full_data = colMeans(holzinger_swineford[names(d)[1:3]]))
#'
#' @export
cov_fiml <- function(data, aux = NULL) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    .dmar_refuse(c(
      "cov_fiml() requires the 'lavaan' package, whose EM algorithm ",
      "computes the estimates. Install it with ",
      "install.packages(\"lavaan\")."), "package_required", package = "lavaan")
  }
  if (!is.data.frame(data) && !is.matrix(data)) {
    .dmar_refuse_argument(
      "'data' must be a data frame or a numeric matrix.", arg = "data")
  }
  # An unnamed matrix gets the column names V1, V2, ... here.
  data <- as.data.frame(data)
  if (nrow(data) == 0L || ncol(data) == 0L) {
    .dmar_refuse_argument("'data' has no rows or no columns.", arg = "data")
  }
  if (anyDuplicated(names(data))) {
    .dmar_refuse_argument(c(
      "'data' has more than one column named ",
      paste(sQuote(unique(names(data)[duplicated(names(data))])),
            collapse = ", "), "."), arg = "data")
  }
  not_numeric <- names(data)[!vapply(data, is.numeric, logical(1L))]
  if (length(not_numeric)) {
    .dmar_refuse_argument(c(
      "cov_fiml() estimates the covariance matrix of continuous ",
      "variables, and these columns of 'data' are not numeric: ",
      paste(sQuote(not_numeric), collapse = ", "), ". Recode a ",
      "variable that is continuous as numeric, and leave out a ",
      "categorical one (the matrix for ordered categorical items is a ",
      "polychoric one, which lavaan::lavCor() estimates with ",
      "'ordered')."), arg = "data")
  }

  if (!is.null(aux)) {
    if (!is.character(aux) || length(aux) < 1L || anyNA(aux)) {
      .dmar_refuse_argument(
        "'aux' must be a character vector of column names in 'data'.",
        arg = c("aux", "data"))
    }
    if (anyDuplicated(aux)) {
      .dmar_refuse_argument("'aux' contains duplicated names.", arg = "aux")
    }
    absent <- setdiff(aux, names(data))
    if (length(absent)) {
      .dmar_refuse_argument(c(
        "Auxiliary variable(s) not found in 'data': ",
        paste(sQuote(absent), collapse = ", "), ". Supply 'data' ",
        "containing both the analysis variables and the auxiliary ",
        "columns; the columns 'aux' does not name are the analysis ",
        "variables."), arg = c("data", "aux"))
    }
  }
  vars <- setdiff(names(data), aux)
  if (length(vars) == 0L) {
    .dmar_refuse_argument(c(
      "Every column of 'data' is named in 'aux', so no analysis ",
      "variable remains; an auxiliary variable informs the estimates ",
      "but is not one of them."), arg = c("data", "aux"))
  }

  cols <- c(vars, aux)
  x <- matrix(vapply(data[cols], as.numeric, numeric(nrow(data))),
              nrow = nrow(data), dimnames = list(NULL, cols))
  if (any(is.infinite(x))) {
    .dmar_refuse_argument(
      "'data' contains infinite values; code a missing value as NA.",
      arg = "data")
  }
  observed <- !is.na(x)
  # A variable needs two distinct observed values for its variance to be
  # estimable; a constant has none, and the matrix would be singular.
  n_distinct <- vapply(seq_len(ncol(x)), function(j) {
    length(unique(x[observed[, j], j]))
  }, integer(1L))
  flat <- cols[n_distinct < 2L]
  if (length(flat)) {
    .dmar_refuse_argument(c(
      "These columns have fewer than two distinct observed values, so ",
      "their variance is zero or not estimable: ",
      paste(sQuote(flat), collapse = ", "), "."), arg = NULL)
  }

  # A case with every variable missing carries no information; lavaan
  # would set it aside with a warning, so it is set aside here first.
  used <- rowSums(observed) > 0L
  x <- x[used, , drop = FALSE]
  observed <- observed[used, , drop = FALSE]
  N <- nrow(x)
  if (N < 2L) {
    .dmar_refuse_argument(
      "Fewer than two cases have an observed value.", arg = NULL)
  }

  # A pair never observed together has a covariance the data do not
  # identify (lavaan only warns, and its EM returns an arbitrary value),
  # and a pair observed together only where one of them takes a single
  # value has one they do not inform (.stop_single_valued_together()).
  .stop_never_observed_together(x)
  .stop_single_valued_together(x)
  .stop_fiml_unbounded(x, "cov_fiml()")

  # lavaan names the variables of a warning from its EM V1, V2, ... (and
  # would name them v1, v2, ..., the names .cov_fiml_em() gives it); the
  # warning is raised again with the names of 'data'.
  em <- withCallingHandlers(.cov_fiml_em(x), warning = function(w) {
    msg <- conditionMessage(w)
    named <- .cov_fiml_user_names(msg, cols)
    if (!identical(named, msg)) {
      .dmar_warn(named, "result_qualified")
      invokeRestart("muffleWarning")
    }
  })
  if (!em$converged) {
    .dmar_warn(c(
      "The EM algorithm reached its iteration limit before meeting ",
      "its convergence criterion, so 'S' and 'mu' are its last ",
      "iterate rather than the maximum likelihood estimates."),
      "result_qualified")
  }

  key <- apply(observed, 1L, function(r) paste(as.integer(r), collapse = ""))
  list(
    S          = em$S[vars, vars, drop = FALSE],
    mu         = em$mu[vars],
    N          = N,
    N_complete = sum(rowSums(observed[, vars, drop = FALSE]) ==
                       length(vars)),
    patterns   = length(unique(key)),
    converged  = em$converged
  )
}


# .cov_fiml_em(x, em_args = NULL)
#
# The EM estimates of the saturated multivariate normal model, computed by
# lavaan (lavaan::lavCor() with missing = "ml" estimates nothing beyond
# the unrestricted h1 moments, which it computes by EM, and output = "fit"
# returns the object that holds them). The columns are renamed v1, v2, ...
# for lavaan, so a column name that is not valid lavaan syntax cannot
# break the call, and the names are restored on the way out.
#
# The EM runs on the variables centered at their available case means and
# divided by their available case SDs, and the estimates are transformed
# back, S = D S_z D and mu = m + D mu_z. lavaan's EM stops when the
# largest change in any mean or covariance is below an absolute 1e-5, and
# it reads an eigenvalue of the covariance matrix below 1e-5 as a sign
# that the matrix is singular, so in the raw metric both depend on the
# units, and lavCor() does not rescale (lavaan's rescale.data applies only
# to a model it fits). Through the first version of this function, which
# ran the EM in the raw metric, variables with SDs near .001 had their EM
# stop after one accelerated step, with correlations off by up to .16,
# converged = TRUE, and a warning that the variances were near zero; the
# error was smaller but still there at SDs of .01 to .1 when much
# information was missing. The ML estimates are equivariant under the
# affine map, so the transformed estimates are the ML estimates of the
# original variables, and the convergence criterion is now relative to
# each variable's scale: a variable multiplied by a constant gives the
# same standardized values, to rounding error, and so the same EM
# iterations, and its estimates change by that constant as the ML
# estimates do.
#
# Also the engine of the reliability family's FIML covariance matrix,
# .fiml_cov() in R/reliability_internals.R.
#
# Inputs:
#   x       : numeric matrix (or data frame of numeric columns) with column
#             names; every row has at least one observed value (cov_fiml()
#             and .fiml_cov() set aside the others), and every pair of
#             columns should be observed together at least once (cov_fiml()
#             refuses a pair that is not; lavaan warns about such a pair,
#             and the EM value of its covariance is then not determined by
#             the data).
#   em_args : NULL for lavaan's defaults, or a list passed as em.h1.args
#             (the tests use it to force the iteration limit).
#
# Returns:
#   A list with S (the covariance matrix, divisor N), mu (the named mean
#   vector), and converged (FALSE when lavaan reported that the EM reached
#   its iteration limit; that warning is muffled here, and cov_fiml()
#   raises its own).
.cov_fiml_em <- function(x, em_args = NULL) {
  original <- colnames(x)
  internal <- paste0("v", seq_along(original))
  x <- matrix(vapply(as.data.frame(x), as.numeric, numeric(NROW(x))),
              nrow = NROW(x))
  center <- colMeans(x, na.rm = TRUE)
  center[!is.finite(center)] <- 0
  scale <- apply(x, 2L, stats::sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale <= 0] <- 1
  z <- sweep(sweep(x, 2L, center), 2L, scale, "/")
  d <- as.data.frame(z)
  names(d) <- internal
  args <- list(d, missing = "ml", output = "fit")
  if (!is.null(em_args)) args$em.h1.args <- em_args

  state <- new.env(parent = emptyenv())
  state$capped <- FALSE
  # lavaan's notice that the EM stopped at its iteration limit is replaced
  # by cov_fiml()'s own warning; every other warning passes through.
  at_limit <- function(w) {
    msg <- conditionMessage(w)
    if (grepl("Maximum number of iterations", msg, fixed = TRUE) &&
        grepl("EM", msg, fixed = TRUE)) {
      state$capped <- TRUE
      invokeRestart("muffleWarning")
    }
  }
  fit <- withCallingHandlers(do.call(lavaan::lavCor, args),
                             warning = at_limit)
  h1 <- lavaan::lavInspect(fit, "sampstat.h1")
  S <- matrix(as.numeric(h1$cov), nrow = length(internal),
              dimnames = dimnames(h1$cov))[internal, internal, drop = FALSE]
  mu <- stats::setNames(as.numeric(h1$mean), names(h1$mean))[internal]
  S <- S * outer(scale, scale)
  mu <- center + scale * mu
  dimnames(S) <- list(original, original)
  names(mu) <- original
  list(S = S, mu = mu, converged = !state$capped)
}


# .stop_never_observed_together(data,
#                               method = "full information maximum likelihood")
#
# Refuses, by name, the pairs of variables in `data` that no case
# observes together. The data carry no information about the covariance
# of such a pair: the likelihood is flat in that direction, lavaan only
# warns that the covariance is not identified, and its EM returns an
# arbitrary value for it, on which anything computed from the matrix then
# rests (in the reliability family, a coefficient with a standard error
# of zero or none, since the information matrix is singular). One
# refusal, in one set of words, serves every caller: cov_fiml() before
# its EM, and the reliability family before a fit by full information
# maximum likelihood or two-stage estimation (.omega_fit_cfa(),
# .fiml_cov(), and the two-stage bootstrap replication
# .omega_fit_cfa_ts_estimates() in R/reliability_internals.R), where it
# runs on the data and on every bootstrap resample and BCa jackknife
# deletion, since thinning can leave a pair with no case in common; there
# the refusal is an error the bootstrap engine counts and quotes. Until
# it moved here, cov_fiml() and the reliability family each carried a
# copy, worded differently.
#
# Inputs:
#   data   : data frame or matrix of the variables an estimate would use
#            (the analysis variables or items, and any auxiliary
#            variables), with column names.
#   method : the words for the method that needs the covariances, as the
#            message names it ("two-stage maximum likelihood" from the
#            reliability family under missing = "two_stage";
#            .relia_incomplete_label()).
#   remedy : the sentence that ends the message in place of the default
#            (leave out one variable of each pair), for a caller whose
#            remedy differs (the fraction of missing information,
#            .fmi_core() in R/missing_information_internals.R).
#   class  : classes the error carries before "error" (NULL for none).
# Returns:
#   NULL, invisibly, when every pair is observed together at least once;
#   otherwise an error.
.stop_never_observed_together <- function(
    data, method = "full information maximum likelihood", remedy = NULL,
    class = NULL) {
  data <- as.data.frame(data)
  coverage <- crossprod(!is.na(data))
  zero <- which(coverage == 0 & upper.tri(coverage), arr.ind = TRUE)
  if (!nrow(zero)) return(invisible(NULL))
  zero <- zero[order(zero[, 1L], zero[, 2L]), , drop = FALSE]
  cols <- colnames(data)
  pairs <- paste0(cols[zero[, 1L]], " with ", cols[zero[, 2L]])
  if (is.null(remedy)) {
    remedy <- paste0(
      "Leave out one variable of each pair (an auxiliary variable where the ",
      "pair has one, since it only informs the estimates); ",
      "missing_coverage() reports the cases behind every pair.")
  }
  .dmar_refuse(
    paste0(
      "These pairs of variables are never observed together, so the ",
      "data carry no information about their covariance, which ",
      method, " needs: ", paste(pairs, collapse = "; "), ". lavaan only ",
      "warns about such a pair, and its estimate of that covariance is ",
      "arbitrary, as is anything computed from it. ", remedy),
    sub("^dmar_", "", class))
}


# .single_valued_pairs(data)
#
# The pairs of variables in `data` that are observed together only in
# rows where one of them takes a single value (which includes a pair one
# row alone observes together). Such rows carry no information about how
# the two vary together. .stop_never_observed_together() lets the pair
# pass, since some row observes both, but the regression of one on the
# other in the rows that observe both is not determined, so the estimate
# of their covariance under the saturated normal model rests on the form
# of the model (the link it makes between a conditional mean and a
# conditional variance) rather than on the data, and when one of the two
# is observed nowhere else it is not determined at all. The case that
# prompted the check: on rows 1 to 100 of four Holzinger and Swineford
# tests with t2_cubes set to 25 in the 90 rows that observe
# t1_visual_perception, and left as it was in the 10 that miss it, full
# information maximum likelihood returned robust omega 1.11 (a Heywood
# case, with t1's error variance at -48 against an observed variance of
# 50) and alpha with no standard error, while listwise deletion, which
# refuses t2 as constant in the complete rows, had offered FIML as the
# remedy.
#
# For each ordered pair (j, k) the question is whether j takes a single
# value in the rows that observe both. The sum and the sum of squares of
# each variable over the rows that observe each other variable come from
# three cross products of the data (centered at the observed means, so
# the sums of squares keep their precision); a pair whose sum of squared
# deviations there is zero to within rounding is a candidate, and each
# candidate is confirmed on its values, so the verdict is exact and the
# cost is that of the cross products, O(n p^2), which matters because the
# check runs on every bootstrap replication.
#
# Inputs:
#   data : data frame or matrix with column names.
# Returns:
#   A matrix with one row per unordered pair and columns j, k (the column
#   indices, j < k), flat_j, flat_k (1 when that variable takes a single
#   value in the rows observing both), and together (the number of those
#   rows); zero rows when there is none.
.single_valued_pairs <- function(data) {
  data <- as.data.frame(data)
  n <- nrow(data)
  p <- ncol(data)
  none <- matrix(numeric(0), 0L, 5L,
                 dimnames = list(NULL, c("j", "k", "flat_j", "flat_k",
                                         "together")))
  if (p < 2L || n == 0L) return(none)
  x <- matrix(vapply(data, as.numeric, numeric(n)), nrow = n)
  obs <- !is.na(x)
  o <- obs * 1
  together <- crossprod(o)
  m <- colMeans(x, na.rm = TRUE)
  m[!is.finite(m)] <- 0
  z <- sweep(x, 2L, m)
  z[!obs] <- 0
  # s1[j, k] and s2[j, k]: the sum and the sum of squares of variable j
  # (centered) over the rows that observe both j and k.
  s1 <- crossprod(z, o)
  s2 <- crossprod(z * z, o)
  ss <- s2 - s1^2 / pmax(together, 1)
  cand <- together > 0 & ss <= 1e-8 * s2
  diag(cand) <- FALSE
  if (!any(cand)) return(none)
  idx <- which(cand, arr.ind = TRUE)
  flat <- vapply(seq_len(nrow(idx)), function(r) {
    v <- x[obs[, idx[r, 1L]] & obs[, idx[r, 2L]], idx[r, 1L]]
    all(v == v[1L])
  }, logical(1L))
  idx <- idx[flat, , drop = FALSE]
  if (!nrow(idx)) return(none)
  j <- pmin(idx[, 1L], idx[, 2L])
  k <- pmax(idx[, 1L], idx[, 2L])
  key <- unique(cbind(j, k))
  key <- key[order(key[, 1L], key[, 2L]), , drop = FALSE]
  flat_of <- function(a, b) {
    as.numeric(vapply(seq_along(a), function(i) {
      any(idx[, 1L] == a[i] & idx[, 2L] == b[i])
    }, logical(1L)))
  }
  out <- cbind(j = key[, 1L], k = key[, 2L],
               flat_j = flat_of(key[, 1L], key[, 2L]),
               flat_k = flat_of(key[, 2L], key[, 1L]),
               together = together[key])
  out
}


# .single_valued_pair_phrase(pairs, cols, r)
#
# Row r of .single_valued_pairs() in words: "t1 with t2 (t2 takes one
# value in the 90 rows that observe both)", or "y1 with y2 (one row
# observes both)".
.single_valued_pair_phrase <- function(pairs, cols, r) {
  j <- pairs[r, "j"]
  k <- pairs[r, "k"]
  n <- pairs[r, "together"]
  paste0(cols[j], " with ", cols[k], " (", if (n == 1) {
    "one row observes both"
  } else if (pairs[r, "flat_j"] == 1 && pairs[r, "flat_k"] == 1) {
    paste0("each takes one value in the ", n, " rows that observe both")
  } else {
    paste0(if (pairs[r, "flat_j"] == 1) cols[j] else cols[k],
           " takes one value in the ", n, " rows that observe both")
  }, ")")
}


# .stop_single_valued_together(data,
#                              method = "full information maximum likelihood")
#
# Refuses, by name, the pairs .single_valued_pairs() finds, as
# .stop_never_observed_together() refuses a pair that no row observes:
# every caller of that function calls this one after it (cov_fiml()
# before its EM; .omega_fit_cfa(), .fiml_cov(), and
# .omega_fit_cfa_ts_estimates() in R/reliability_internals.R before a fit
# by full information maximum likelihood or two-stage estimation, on the
# data and on every bootstrap resample and BCa jackknife deletion, where
# the refusal is an error the bootstrap engine counts and quotes).
#
# Inputs:
#   data, method : as in .stop_never_observed_together().
# Returns:
#   NULL, invisibly, when there is no such pair; otherwise an error.
.stop_single_valued_together <- function(
    data, method = "full information maximum likelihood") {
  pairs <- .single_valued_pairs(data)
  if (!nrow(pairs)) return(invisible(NULL))
  cols <- colnames(as.data.frame(data))
  described <- vapply(seq_len(nrow(pairs)), function(r) {
    .single_valued_pair_phrase(pairs, cols, r)
  }, character(1L))
  .dmar_refuse_argument(c(
    "These pairs of variables are observed together only in rows ",
    "where one of them takes a single value, so those rows carry no ",
    "information about how the two vary together, which ", method,
    " needs to estimate their covariance: ",
    paste(described, collapse = "; "), ". Its estimate would rest on ",
    "the form of the normal model rather than on the data, or not be ",
    "determined at all, and so would anything computed from it. Leave ",
    "out one variable of each pair (an auxiliary variable where the ",
    "pair has one, since it only informs the estimates); ",
    "missing_coverage() reports the cases behind every pair."), arg = NULL)
}


# .cov_fiml_user_names(msg, names)
#
# A lavaan message with the variables it names given the names they have
# in 'data'. lavaan calls the columns of a matrix without names V1, V2,
# ... (as in the warning its EM raises when the covariance matrix is near
# singular), and .cov_fiml_em() hands it the columns as v1, v2, ...; a
# token of either form whose number is a column of `names` is replaced,
# and anything else is left as it is.
.cov_fiml_user_names <- function(msg, names) {
  hits <- gregexpr("\\b[Vv][0-9]+\\b", msg, perl = TRUE)
  regmatches(msg, hits) <- lapply(regmatches(msg, hits), function(tok) {
    k <- suppressWarnings(as.integer(substring(tok, 2L)))
    ok <- !is.na(k) & k >= 1L & k <= length(names)
    tok[ok] <- names[k[ok]]
    tok
  })
  msg
}
