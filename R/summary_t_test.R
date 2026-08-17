# Two-sample t test from summary statistics.
#' Two-Sample \emph{t} Test From Summary Statistics
#'
#' Computes a two-sample \emph{t} test (pooled or Welch) directly from
#' the per-group means, standard deviations, and sample sizes, without
#' requiring access to the raw observations. Returns the test statistic,
#' degrees of freedom, \emph{p}-value, and a CI on the mean difference
#' in a \code{data.frame}. Useful for re-analyses from published
#' papers that report only the summary numbers.
#'
#' @param mean_1,mean_2 Group sample means.
#' @param sd_1,sd_2 Group sample standard deviations.
#' @param n_1,n_2 Group sample sizes.
#' @param mu Null value of the mean difference \eqn{\mu_1 - \mu_2}.
#'   Default \code{0}.
#' @param var_equal Logical. If \code{TRUE} (default), uses Student's
#'   pooled-variance \emph{t}. If \code{FALSE}, uses Welch's separate-
#'   variance \emph{t} with Satterthwaite degrees of freedom.
#' @param alternative One of \code{"two_sided"} (default; the base-R
#'   spelling \code{"two.sided"} is accepted as an alias),
#'   \code{"less"}, or \code{"greater"}.
#' @param conf_level Confidence level for the CI on the mean
#'   difference. Default \code{0.95}.
#'
#' @return A \code{data.frame} with rows for the mean difference,
#'   the \emph{t} statistic, degrees of freedom, \emph{p}-value, and
#'   the CI lower and upper limits on the mean difference.
#'
#' @details
#' \strong{Pooled-variance \emph{t} (Student, 1908).} Under
#' \eqn{\sigma_1 = \sigma_2}, the pooled SD is
#' \eqn{s_p = \sqrt{((n_1 - 1) s_1^2 + (n_2 - 1) s_2^2) / (n_1 + n_2 - 2)}},
#' the test statistic is
#' \eqn{t = (\bar x_1 - \bar x_2 - \mu_0) / (s_p \sqrt{1 / n_1 + 1 / n_2})},
#' and \eqn{df = n_1 + n_2 - 2}.
#'
#' \strong{Welch's \emph{t} (Welch, 1947).} Under unequal variances,
#' \eqn{t = (\bar x_1 - \bar x_2 - \mu_0) /
#'   \sqrt{s_1^2 / n_1 + s_2^2 / n_2}}
#' with Satterthwaite degrees of freedom (see \code{\link{welch_t}}).
#'
#' \strong{Choosing pooled vs Welch.} Methodological reviews now
#' recommend Welch as the default (Delacre, Lakens, & Leys, 2017;
#' Ruxton, 2006). Pooled-variance \emph{t} is preserved here primarily
#' for reproducing analyses from older sources that used it.
#'
#' @references
#' Delacre, M., Lakens, D., & Leys, C. (2017). Why psychologists should
#'   by default use Welch's \emph{t}-test instead of Student's
#'   \emph{t}-test. \emph{International Review of Social Psychology,
#'   30}(1), 92--101. \doi{10.5334/irsp.82}
#'
#' Ruxton, G. D. (2006). The unequal variance \emph{t}-test is an
#'   underused alternative to Student's \emph{t}-test and the
#'   Mann-Whitney \emph{U} test. \emph{Behavioral Ecology, 17}(4),
#'   688--690. \doi{10.1093/beheco/ark016}
#'
#' Snedecor, G. W., & Cochran, W. G. (1989). \emph{Statistical methods}
#'   (8th ed.). Iowa State University Press.
#'
#' Student. (1908). The probable error of a mean. \emph{Biometrika, 6}(1),
#'   1--25. \doi{10.2307/2331554}
#'
#' Welch, B. L. (1947). The generalization of "Student's" problem when
#'   several different population variances are involved.
#'   \emph{Biometrika, 34}(1/2), 28--35.
#'
#' @seealso \code{\link{welch_t}}, \code{\link[stats]{t.test}},
#'   \code{\link{smd}}
#'
#' @examples
#' # 1. Re-analysis from published summary statistics:
#' #        Group A: M = 100, SD = 15, n = 30
#' #        Group B: M = 108, SD = 18, n = 25
#' summary_t_test(mean_1 = 100, sd_1 = 15, n_1 = 30,
#'                mean_2 = 108, sd_2 = 18, n_2 = 25)
#'
#' # 2. Welch version for the same data:
#' summary_t_test(mean_1 = 100, sd_1 = 15, n_1 = 30,
#'                mean_2 = 108, sd_2 = 18, n_2 = 25,
#'                var_equal = FALSE)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family hypothesis tests
#'
#' @export

summary_t_test <- function(mean_1, sd_1, n_1,
                           mean_2, sd_2, n_2,
                           mu = 0,
                           var_equal = TRUE,
                           alternative = c("two_sided", "less", "greater"),
                           conf_level = 0.95) {
  alternative <- .match_alternative(alternative)

  for (nm in c("mean_1", "mean_2", "sd_1", "sd_2", "n_1", "n_2", "mu")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L)
      stop(sprintf("'%s' must be a single numeric value.", nm))
  }
  if (sd_1 < 0 || sd_2 < 0)
    stop("'sd_1' and 'sd_2' must be non-negative.")
  if (n_1 < 2 || n_2 < 2)
    stop("'n_1' and 'n_2' must each be at least 2.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  diff <- mean_1 - mean_2

  if (var_equal) {
    s2_p <- ((n_1 - 1) * sd_1^2 + (n_2 - 1) * sd_2^2) / (n_1 + n_2 - 2)
    se   <- sqrt(s2_p * (1 / n_1 + 1 / n_2))
    df   <- n_1 + n_2 - 2
  } else {
    se <- sqrt(sd_1^2 / n_1 + sd_2^2 / n_2)
    df <- (sd_1^2 / n_1 + sd_2^2 / n_2)^2 /
          ((sd_1^2 / n_1)^2 / (n_1 - 1) +
           (sd_2^2 / n_2)^2 / (n_2 - 1))
  }
  if (se == 0)
    stop("Standard error of the mean difference is zero; t is undefined.")

  t_stat <- (diff - mu) / se
  p_value <- switch(alternative,
    two_sided = 2 * stats::pt(-abs(t_stat), df),
    less      = stats::pt(t_stat, df),
    greater   = stats::pt(t_stat, df, lower.tail = FALSE)
  )

  alpha <- 1 - conf_level
  t_crit <- switch(alternative,
    two_sided = stats::qt(1 - alpha / 2, df),
    less      = stats::qt(conf_level,    df),
    greater   = stats::qt(conf_level,    df)
  )
  ci_lo <- switch(alternative,
    two_sided = diff - t_crit * se,
    less      = -Inf,
    greater   = diff - t_crit * se
  )
  ci_hi <- switch(alternative,
    two_sided = diff + t_crit * se,
    less      = diff + t_crit * se,
    greater   = Inf
  )

  out <- data.frame(
    term  = c("mean_difference", "t_statistic", "df", "p_value",
              "lower_limit", "upper_limit"),
    value = c(diff, t_stat, df, p_value, ci_lo, ci_hi),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out <- .as_dmar_tbl(out, conf_level = conf_level, p_terms = "p_value")
  # Leading subclass for broom dispatch, layered over the dmar_tbl class so
  # print() still falls through to print.dmar_tbl (see R/dmar_tidiers.R).
  class(out) <- c("dmar_summary_t_test", class(out))
  out
}


#' Broom-Style Tidy / Glance Methods for \code{summary_t_test()}
#'
#' \code{tidy()} returns the single mean-difference estimate and its
#' confidence interval in the \pkg{broom} convention; \code{glance()}
#' coincides with it, since a two-sample \emph{t} test reports one
#' estimand and there are no extra model-level statistics to add.
#'
#' @param x A \code{dmar_summary_t_test} object returned by
#'   \code{\link{summary_t_test}}.
#' @param \dots Unused.
#'
#' @return A one-row \code{data.frame} with columns \code{term},
#'   \code{estimate}, \code{ci_lower}, \code{ci_upper}, \code{statistic},
#'   \code{df}, \code{p_value}, and \code{conf_level}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @name summary_t_test_broom
NULL

#' @rdname summary_t_test_broom
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_summary_t_test <- function(x, ...) {
  pick <- function(term) x$value[x$term == term]
  data.frame(
    term       = "mean_difference",
    estimate   = pick("mean_difference"),
    ci_lower   = pick("lower_limit"),
    ci_upper  = pick("upper_limit"),
    statistic  = pick("t_statistic"),
    df         = pick("df"),
    p_value    = pick("p_value"),
    conf_level = attr(x, "conf_level") %||% NA_real_,
    stringsAsFactors = FALSE,
    row.names  = NULL
  )
}

#' @rdname summary_t_test_broom
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_summary_t_test <- function(x, ...) {
  tidy.dmar_summary_t_test(x, ...)
}
