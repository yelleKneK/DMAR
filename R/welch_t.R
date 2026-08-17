# Welch's separate-variance t test in tidy form.
#' Welch's Separate-Variance \emph{t} Test
#'
#' Computes Welch's (1947) separate-variance \emph{t} test, the
#' Satterthwaite (1946) approximation for unequal-variance two-sample
#' inference, and returns the test statistic, Satterthwaite degrees of
#' freedom, \emph{p}-value, point estimate of the mean difference, and
#' a confidence interval on the mean difference, all in a tidy
#' \code{data.frame}. Unlike Student's pooled-variance \emph{t} test
#' (\code{stats::t.test(..., var.equal = TRUE)}), Welch's test does
#' \emph{not} assume the two populations have equal variances, and
#' should be the default choice in applied work (Delacre, Lakens, &
#' Leys, 2017).
#'
#' @param x,y Numeric vectors of observations from the two groups.
#'   The two groups are independent and need not be the same length;
#'   \code{NA}s are removed from each vector separately.
#' @param mu Null value of the mean difference \eqn{\mu_1 - \mu_2}.
#'   Default \code{0}.
#' @param alternative One of \code{"two_sided"} (default; the base-R
#'   spelling \code{"two.sided"} is accepted as an alias),
#'   \code{"less"}, or \code{"greater"}, defining the direction of the
#'   alternative hypothesis.
#' @param conf_level Confidence level for the CI on the mean
#'   difference. Default \code{0.95}.
#'
#' @return A \code{data.frame} with rows for the mean difference
#'   \eqn{\bar x - \bar y}, the Welch \emph{t}-statistic, Satterthwaite
#'   degrees of freedom, \emph{p}-value, the CI lower and upper limits
#'   on the mean difference, and the per-group means, SDs, and \emph{n}.
#'
#' @details
#' \strong{Test statistic.} Welch's \emph{t} is
#' \deqn{t \;=\; \frac{\bar x - \bar y - \mu_0}
#'                    {\sqrt{s_1^2 / n_1 + s_2^2 / n_2}},}
#' which is referred to a \emph{t} distribution on the Satterthwaite
#' (1946) approximate degrees of freedom
#' \deqn{df \;=\; \frac{(s_1^2 / n_1 + s_2^2 / n_2)^2}
#'                     {(s_1^2 / n_1)^2 / (n_1 - 1) +
#'                      (s_2^2 / n_2)^2 / (n_2 - 1)}.}
#'
#' \strong{Why Welch by default.} Student's pooled-variance \emph{t}
#' assumes \eqn{\sigma_1 = \sigma_2}; when that assumption fails it has
#' both inflated and deflated Type I error rates depending on the
#' \eqn{n_1 : n_2} ratio (Ruxton, 2006). Welch's test maintains nominal
#' Type I error across virtually all combinations of \eqn{\sigma_1 /
#' \sigma_2} and \eqn{n_1 / n_2}, with no meaningful loss of power
#' when variances are equal. The American Statistical Association and
#' multiple methodological reviews now recommend Welch as the default
#' (Delacre et al., 2017; Lakens, 2015).
#'
#' \strong{Relation to \code{stats::t.test()}.} The numerical results
#' here match \code{stats::t.test(x, y, var.equal = FALSE)} to machine
#' precision; this function differs only in returning a tidy
#' \code{data.frame} that composes with the rest of DMAR.
#'
#' @references
#' Delacre, M., Lakens, D., & Leys, C. (2017). Why psychologists should
#'   by default use Welch's \emph{t}-test instead of Student's
#'   \emph{t}-test. \emph{International Review of Social Psychology,
#'   30}(1), 92--101. \doi{10.5334/irsp.82}
#'
#' Lakens, D. (2015, January). Always use Welch's t-test instead of
#'   Student's t-test [Blog post]. The 20\% Statistician.
#'   \url{https://daniellakens.blogspot.com/2015/01/always-use-welchs-t-test-instead-of.html}
#'
#' Ruxton, G. D. (2006). The unequal variance \emph{t}-test is an
#'   underused alternative to Student's \emph{t}-test and the
#'   Mann-Whitney \emph{U} test. \emph{Behavioral Ecology, 17}(4),
#'   688--690. \doi{10.1093/beheco/ark016}
#'
#' Satterthwaite, F. E. (1946). An approximate distribution of
#'   estimates of variance components. \emph{Biometrics Bulletin, 2}(6),
#'   110--114.
#'
#' Welch, B. L. (1947). The generalization of "Student's" problem when
#'   several different population variances are involved.
#'   \emph{Biometrika, 34}(1/2), 28--35.
#'
#' @seealso \code{\link[stats]{t.test}}, \code{\link{summary_t_test}},
#'   \code{\link{smd}}, \code{\link{ci_smd}}
#'
#' @examples
#' # 1. Two groups with different variances:
#' set.seed(113)
#' x <- rnorm(20, mean = 100, sd = 15)
#' y <- rnorm(20, mean = 110, sd = 25)
#' welch_t(x, y)
#'
#' # 2. One-sided test:
#' welch_t(x, y, alternative = "less")
#'
#' # 3. Side-by-side comparison with base R's stats::t.test().
#' # The two functions implement the same Welch / Satterthwaite test, so
#' # the t-statistic, Satterthwaite degrees of freedom, p-value, and the
#' # CI on the mean difference match exactly. welch_t() differs only in
#' # what it returns: a data.frame(term, value) rather than a list-
#' # like htest object. The return composes with dplyr / ggplot2
#' # pipelines and avoids stringly-typed access like $statistic.
#' set.seed(113)
#' a <- rnorm(15, mean = 0,   sd = 1)
#' b <- rnorm(20, mean = 0.5, sd = 2)
#'
#' # DMAR (data.frame):
#' dmar_res <- welch_t(a, b, conf_level = 0.95)
#' dmar_res
#'
#' # Base R (htest list):
#' base_res <- stats::t.test(a, b, var.equal = FALSE, conf.level = 0.95)
#' base_res
#'
#' # Verify the four key statistics agree numerically:
#' pick <- function(term) dmar_res$value[dmar_res$term == term]
#' stopifnot(
#'   all.equal(pick("t_statistic"), unname(base_res$statistic)),
#'   all.equal(pick("df"),          unname(base_res$parameter)),
#'   all.equal(pick("p_value"),     base_res$p.value),
#'   all.equal(pick("lower_limit"), base_res$conf.int[1]),
#'   all.equal(pick("upper_limit"), base_res$conf.int[2])
#' )
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family hypothesis tests
#'
#' @export

welch_t <- function(x, y, mu = 0,
                    alternative = c("two_sided", "less", "greater"),
                    conf_level = 0.95) {
  alternative <- .match_alternative(alternative)

  if (!is.numeric(x) || !is.numeric(y))
    stop("'x' and 'y' must be numeric vectors.")
  if (!is.numeric(mu) || length(mu) != 1L)
    stop("'mu' must be a single numeric value.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  n_1 <- length(x); n_2 <- length(y)
  if (n_1 < 2L || n_2 < 2L)
    stop("Each group needs at least 2 non-missing observations.")

  m_1 <- mean(x);   m_2 <- mean(y)
  v_1 <- stats::var(x); v_2 <- stats::var(y)
  if (v_1 == 0 && v_2 == 0)
    stop("Both groups have zero within-group variance; Welch's t is undefined.")

  se   <- sqrt(v_1 / n_1 + v_2 / n_2)
  diff <- m_1 - m_2
  t_stat <- (diff - mu) / se
  df <- (v_1 / n_1 + v_2 / n_2)^2 /
        ((v_1 / n_1)^2 / (n_1 - 1) + (v_2 / n_2)^2 / (n_2 - 1))

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
              "lower_limit", "upper_limit",
              "mean_x", "mean_y", "sd_x", "sd_y", "n_x", "n_y"),
    value = c(diff, t_stat, df, p_value,
              ci_lo, ci_hi,
              m_1, m_2, sqrt(v_1), sqrt(v_2), n_1, n_2),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out <- .as_dmar_tbl(out, conf_level = conf_level, p_terms = "p_value")
  # Leading subclass for broom dispatch, layered over the dmar_tbl class so
  # print() still falls through to print.dmar_tbl (see R/dmar_tidiers.R).
  class(out) <- c("dmar_welch_t", class(out))
  out
}


#' Broom-Style Tidy / Glance Methods for \code{welch_t()}
#'
#' \code{tidy()} returns the single mean-difference estimate and its
#' confidence interval in the \pkg{broom} convention; \code{glance()}
#' coincides with it, since a two-sample \emph{t} test reports one
#' estimand and there are no extra model-level statistics to add.
#'
#' @param x A \code{dmar_welch_t} object returned by \code{\link{welch_t}}.
#' @param \dots Unused.
#'
#' @return A one-row \code{data.frame} with columns \code{term},
#'   \code{estimate}, \code{ci_lower}, \code{ci_upper}, \code{statistic},
#'   \code{df}, \code{p_value}, and \code{conf_level}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @name welch_t_broom
NULL

#' @rdname welch_t_broom
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_welch_t <- function(x, ...) {
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

#' @rdname welch_t_broom
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_welch_t <- function(x, ...) {
  tidy.dmar_welch_t(x, ...)
}
