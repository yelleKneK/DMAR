#' Kuder-Richardson Formula 20 (KR-20) With a Confidence Interval
#'
#' @description
#' Estimates Kuder-Richardson formula 20 (Kuder & Richardson, 1937) for a
#' homogeneous composite scored on dichotomous (0/1) items, and returns a
#' confidence interval for the population coefficient.
#'
#' @details
#' Kuder and Richardson's (1937) formula 20 for a \eqn{J}-item composite
#' of binary items is
#' \deqn{KR_{20} = \frac{J}{J-1}\left(1 - \frac{\sum_{j} p_{j} q_{j}}{s_{Y}^{2}}\right),}
#' where \eqn{p_{j}} is the proportion of respondents endorsing item
#' \emph{j} (i.e., scoring 1), \eqn{q_{j} = 1 - p_{j}}, and
#' \eqn{s_{Y}^{2}} is the variance of the composite score. For
#' dichotomous items \eqn{p_{j} q_{j}} is the item variance, so KR-20 is
#' algebraically identical to coefficient \eqn{\alpha} (Guttman, 1945;
#' Cronbach, 1951) computed from the item covariance matrix. KR-20
#' predates coefficient \eqn{\alpha} by 14 years and Feldt's (1965)
#' \emph{F}-distribution interval was derived specifically for the
#' sampling distribution of KR-20.
#'
#' Because KR-20 is a special case of coefficient \eqn{\alpha} (Guttman,
#' 1945; Cronbach, 1951), the same considerations apply: KR-20 equals
#' the population reliability of the composite under essential
#' \eqn{\tau}-equivalence (i.e., equal factor loadings) and serves as a
#' lower bound otherwise. For dichotomous items see also
#' \code{\link{reliability_omega_categorical}} (categorical omega), which treats
#' the items via a probit-link single-factor model and does not assume
#' equal loadings.
#'
#' KR-20, like \eqn{\alpha}, is at most one but can fall below zero, when
#' the item covariances are negative on balance; on dichotomous items
#' that can arise when an item is scored in the opposite direction to
#' the others (a reversed key). A reliability coefficient cannot be
#' negative, so a negative estimate is returned as computed with a
#' warning that says so, and each confidence limit is clamped into
#' [0, 1]; the interval around an estimate far enough below zero is then
#' [0, 0], which says that the data support no positive reliability for
#' the items as keyed.
#'
#' Available confidence interval methods (set via \code{ci_method}) are
#' the same as for \code{\link{reliability_alpha}}; see that function's
#' \emph{Details} for full descriptions. The default for KR-20 is
#' \code{"feldt"}, the \emph{F}-distribution interval originally
#' developed for KR-20.
#'
#' From two complete rows KR-20 is computed with Feldt's interval, which
#' is defined for any \emph{N} of at least 2 (every item must differ
#' between the two rows, since an item they answer alike does not vary,
#' and the two totals must differ, since a total that does not vary has
#' no reliability).
#' Feldt's \emph{F} reference is exact under multivariate normality and
#' parallel items, which 0/1 items cannot meet, and for them it is the
#' approximation Feldt (1965) derived. The intervals that are not defined
#' there (\code{"bonett"}, \code{"fisher"}, \code{"ml"}, \code{"adf"},
#' their \code{_logistic} variants, and the bootstrap) are refused, as
#' \code{\link{reliability_alpha}} explains, and \code{"hakstian_whalen"},
#' which is defined, has an upper limit of one there whatever the data, at
#' any confidence level above about .90.
#'
#' The menu also includes the nonparametric bootstrap intervals
#' \code{"percentile"}, \code{"bca"}, \code{"bootstrap_se"}, and
#' \code{"bootstrap_se_logistic"}, which resample the rows of
#' \code{data} with replacement \code{B} times and recompute KR-20 on
#' each replication (Efron & Tibshirani, 1993). They are worth the cost
#' when the closed forms are least trustworthy, which for dichotomous
#' items means highly unbalanced item difficulties or a sample size too
#' small for the normal-theory derivations behind \code{"feldt"} and
#' \code{"bonett"}. No bootstrap runs unless \code{ci_method} asks for
#' one; when it does, the default is \code{B = 10000} replications, and
#' supplying \code{seed} makes the interval reproducible. Replications on
#' which KR-20 cannot be computed are set aside and counted, among them a
#' replication in which an item takes a single value, refused by name as
#' the data would be: a resample leaves out a given respondent about
#' 37 percent of the time (\eqn{(1 - 1/N)^N}, near \eqn{e^{-1}}), so an
#' item passed by one respondent loses its variation in about a third of
#' the resamples. The
#' interval is computed from the replications that return a value; the
#' counts travel in the \code{bootstrap} attribute described under
#' \code{\link{reliability_alpha}}, print under the table, and are
#' stated in a warning whenever a replication was set aside.
#'
#' \strong{Comparison with other packages.} The \pkg{psych} package
#' computes the same quantity via \code{\link[psych]{alpha}} (since
#' \eqn{\alpha} on 0/1 data \emph{is} KR-20). \code{reliability_kr20}
#' restricts input to raw 0/1 data so the dichotomous-items assumption
#' cannot be quietly violated, presents the historical formula in the
#' documentation, and accompanies the point estimate with a confidence
#' interval drawn from the methods compared in Kelley and Pornprasertmanit
#' (2016).
#'
#' @param data A numeric matrix or data frame of 0/1 item scores (rows
#'   are respondents, columns are items). Rows with any missing values
#'   are listwise-deleted. Non-binary values trigger an error, and so do
#'   fewer than two complete rows, from which no variance can be
#'   computed. A column that is not numeric, and
#'   an item that takes one value in every complete row (one no
#'   respondent passes, or every respondent does), are refused by name:
#'   such an item adds a constant to the total score, whose reliability
#'   it leaves as it is, while KR-20 would count it among the \eqn{J}
#'   items, at \eqn{J(J-2)/(J-1)^2} times the KR-20 of the other items,
#'   the coefficient of the same total. So is a total score that takes
#'   one value in every complete
#'   row (on two rows, as many items passed in one row as in the other),
#'   whose reliability is not defined.
#'   For incomplete data,
#'   \code{\link{reliability_alpha}} with \code{missing = "fiml"}
#'   computes the same coefficient (alpha of 0/1 items is KR-20) from the
#'   incomplete rows as well. SPSS user-defined missing values, which
#'   \code{haven::read_sav(user_na = TRUE)} keeps as values of a
#'   \code{haven_labelled_spss} column, are missing, as SPSS treats
#'   them.
#' @param ci_method Method for constructing the confidence interval; see
#'   \code{\link{reliability_alpha}} for the full list. Defaults to
#'   \code{"feldt"}.
#' @param conf_level Confidence level for the interval (1 - Type I error
#'   rate). Defaults to \code{0.95}.
#' @param B Number of bootstrap replications when a bootstrap method is
#'   selected, a single whole number. Defaults to \code{10000}. The fewest
#'   accepted depends on the interval method and the level: at the default
#'   \code{conf_level}, 400 for \code{"percentile"} (10 replications beyond
#'   each limit), 2000 for \code{"bca"} (50 beyond each nominal limit,
#'   because the adjusted limits sit farther into the tails), and 100 for
#'   \code{"bootstrap_se"} and \code{"bootstrap_se_logistic"} (a standard
#'   deviation of the replications); at 0.99 the first two are 2000 and
#'   10000. A smaller \code{B} is refused before anything is computed
#'   (class \code{dmar_count_invalid}). The count is held to the number
#'   requested.
#'   Replications set aside afterward are counted and reported in a warning
#'   of class \code{dmar_bootstrap_replications_dropped}, which names the
#'   function and says to request more when those that return a value are
#'   fewer than the number required; when fewer than half of the \code{B}
#'   replications return a value, the interval is withheld (limits and
#'   bootstrap standard error \code{NA}) under a warning of class
#'   \code{dmar_bootstrap_insufficient}. See \code{\link{dmar_bootstrap}},
#'   which also describes \code{options(dmar.bootstrap_exploratory = TRUE)}
#'   for a first look with a smaller \code{B}.
#' @param seed Random number seed used for bootstrap reproducibility.
#'   Defaults to \code{NULL}, under which the resamples are drawn from,
#'   and advance, the current random number stream; supply an integer
#'   for a reproducible interval, which is set for the call only, the
#'   caller's random number state being restored when it returns.
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
#'
#' @return A \code{data.frame} with columns \code{term} and \code{value}
#' and rows
#' \code{"estimate"} (sample KR-20),
#' \code{"se"} (standard error on the coefficient scale, \code{NA} for
#' methods that do not produce one; for the transformation-based
#' intervals \code{"fisher"}, \code{"bonett"}, and
#' \code{"hakstian_whalen"} it is the delta method back-transform of the
#' transformation-scale standard error),
#' \code{"se_transformed"} (only for those transformation-based
#' intervals: the standard error on the transformation scale, with the
#' scale named in the attribute \code{se_transform_scale}:
#' \code{"fisher_z"}, \code{"log(1-alpha)"}, or \code{"cube_root"}),
#' \code{"ci_lower"} and \code{"ci_upper"} (each clamped into
#' [0, 1], so the limits stay in order),
#' \code{"conf_level"}, \code{"N"} (effective sample size after
#' listwise deletion), \code{"N_complete"} (the complete cases; equal
#' to \code{"N"} here, and carried so the whole reliability family
#' returns one shape), and \code{"n_items"} (number of items). Attributes
#' \code{coefficient} (\code{"kr20"}) and \code{ci_method} record the
#' computation; bootstrap calls also record \code{B}, the replication
#' counts in the \code{bootstrap} attribute (see
#' \code{\link{reliability_alpha}}; KR-20 fits no model, so no
#' replication is improper), and the number of replications the
#' interval uses in the attribute \code{B_used}, with
#' \code{n_dropped}, \code{B_required}, and \code{bootstrap_status}
#' (see \code{\link{dmar_bootstrap}}).
#'
#' @references
#' Cronbach, L. J. (1951). Coefficient alpha and the internal structure
#'   of tests. \emph{Psychometrika, 16}(3), 297--334.
#'
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Feldt, L. S. (1965). The approximate sampling distribution of
#'   Kuder-Richardson reliability coefficient twenty.
#'   \emph{Psychometrika, 30}, 357--370.
#'
#' Guttman, L. (1945). A basis for analyzing test-retest reliability.
#'   \emph{Psychometrika, 10}(4), 255--282.
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
#' Kuder, G. F., & Richardson, M. W. (1937). The theory of the estimation
#'   of test reliability. \emph{Psychometrika, 2}, 151--160.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
#'   reliability coefficients: Accuracy in parameter estimation via narrow
#'   confidence intervals. \emph{British Journal of Mathematical and
#'   Statistical Psychology, 65}, 371--401.
#'   \doi{10.1111/j.2044-8317.2011.02030.x}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{reliability}} (general wrapper),
#' \code{\link{reliability_alpha}},
#' \code{\link{reliability_omega_categorical}},
#' \code{\link[psych]{alpha}}.
#'
#' @examples
#' set.seed(113)
#' # Ten dichotomous items with a single underlying ability.
#' N <- 300
#' J <- 10
#' ability <- rnorm(N)
#' loadings <- rep(0.6, J)
#' latent <- outer(ability, loadings) +
#'           matrix(rnorm(N * J, sd = sqrt(1 - 0.6^2)), N, J)
#' items <- (latent > 0) * 1
#' colnames(items) <- paste0("y", seq_len(J))
#'
#' reliability_kr20(data = items)
#' reliability_kr20(data = items, ci_method = "bonett")
#'
#' # A bootstrap interval recomputes KR-20 on each of B resamples of the
#' # rows. B = 500 keeps the example quick; a reported interval deserves
#' # the default B = 10000, and seed makes the interval reproducible.
#' reliability_kr20(data = items, ci_method = "percentile", B = 500,
#'                  seed = 113)
#'
#' @keywords htest multivariate
#' @family reliability
#'
#' @export

reliability_kr20 <- function(data,
                             ci_method = c("feldt", "bonett", "fisher",
                                           "hakstian_whalen",
                                           "ml", "ml_logistic",
                                           "adf", "adf_logistic",
                                           "bootstrap_se",
                                           "bootstrap_se_logistic",
                                           "percentile", "bca", "none"),
                             conf_level = 0.95,
                             B = 10000,
                             seed = NULL,
                             ncpus = 1L,
                             progress = NULL) {
  .dmar_check_conf_level(conf_level)
  ci_method <- match.arg(ci_method)
  ncpus <- .check_ncpus(ncpus)
  progress <- .check_progress(progress)
  if (ci_method %in% c("bootstrap_se", "bootstrap_se_logistic",
                       "percentile", "bca")) {
    .check_boot_B(B, ci_method, conf_level)
  }

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  }
  if (missing(data) || is.null(data)) {
    .dmar_refuse_argument(
      "'data' (a 0/1 matrix or data frame of item scores) is required.",
      arg = "data")
  }

  # SPSS user-defined missing values are missing (.relia_declared_missing()).
  data <- .relia_declared_missing(as.data.frame(data))
  # A column that is not numeric (a factor or character coding of 0 and 1
  # passed the check of the values below as text and stopped stats::cov())
  # is refused by name, as in the rest of the family.
  .relia_check_numeric_items(data)
  complete <- stats::complete.cases(data)
  data_mat <- as.matrix(data[complete, , drop = FALSE])
  if (!all(data_mat %in% c(0, 1))) {
    .dmar_refuse_argument(
      "KR-20 requires items scored 0/1; non-binary values were found.",
      arg = NULL)
  }
  if (ncol(data_mat) < 2L) {
    .dmar_refuse_argument("At least two items are required.", arg = NULL)
  }
  # Fewer than two complete rows leave no variance to compute, and the
  # call stops here, with the family's message (offering full information
  # maximum likelihood through reliability_alpha(), since KR-20 has no
  # missing argument). Through version 1.0.1 it returned NA for the
  # estimate and every limit with two warnings that NaNs were produced.
  if (sum(complete) < 2L) {
    .stop_relia_too_few_complete_rows(data, kr20 = TRUE)
  }
  # An item every complete row answers the same way (one nobody passes, or
  # everybody does) is refused by name, as alpha refuses it: KR-20 with it
  # counts it among the J items, J(J - 2)/(J - 1)^2 times KR-20 without
  # it, the coefficient of the same total score (through version 1.0.1 it
  # was counted).
  .relia_check_constant_items(as.data.frame(data_mat), "listwise",
                              all_rows = data, kr20 = TRUE)

  N <- nrow(data_mat)
  J <- ncol(data_mat)
  S <- stats::cov(data_mat)
  estimate <- .kr20_from_data(data_mat)

  # An interval that refuses two rows names Feldt's, which approximates
  # the sampling distribution of KR-20 (.relia_two_rows_alternative()).
  two_rows <- .relia_two_rows_alternative(kr20 = TRUE)
  ci <- switch(
    ci_method,
    none = list(se = NA_real_, lower = NA_real_, upper = NA_real_),
    feldt = .ci_feldt(estimate, N, J, conf_level),
    fisher = .ci_fisher(estimate, N, conf_level),
    bonett = .ci_bonett(estimate, N, J, conf_level, two_rows = two_rows),
    hakstian_whalen = .ci_hakstian_whalen(estimate, N, J, conf_level),
    ml = .ci_alpha_ml(estimate, S, N, J, conf_level, logistic = FALSE,
                      two_rows = two_rows),
    ml_logistic = .ci_alpha_ml(estimate, S, N, J, conf_level,
                               logistic = TRUE, two_rows = two_rows),
    adf = .ci_alpha_adf(estimate, data_mat, conf_level, logistic = FALSE,
                        two_rows = two_rows),
    adf_logistic = .ci_alpha_adf(estimate, data_mat, conf_level,
                                 logistic = TRUE, two_rows = two_rows),
    bootstrap_se = ,
    bootstrap_se_logistic = ,
    percentile = ,
    bca = .bootstrap_ci(
      fun = "reliability_kr20",
      data = data_mat,
      # A replication in which an item takes a single value (a rare item
      # whose passers the resample left out) is refused, as the data are
      # (.relia_checked_point()).
      point_fn = .relia_checked_point(function(d) .kr20_from_data(d)),
      B = B, conf_level = conf_level, kind = ci_method, seed = seed,
      ncpus = ncpus, two_rows = two_rows, progress = progress
    )
  )

  .relia_result(
    estimate = estimate, se = ci$se,
    lower = ci$lower, upper = ci$upper,
    conf_level = conf_level, N = N, J = J,
    coefficient = "kr20",
    ci_method = ci$ci_method %||% ci_method,
    bootstrap = ci$bootstrap,
    B = if (ci_method %in% c("bootstrap_se", "bootstrap_se_logistic",
                             "percentile", "bca")) B else NA_integer_,
    se_transformed = ci$se_transformed,
    se_transform_scale = ci$se_transform_scale
  )
}
