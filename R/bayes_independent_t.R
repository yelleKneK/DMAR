#' Bayesian Independent-Samples \emph{t} Analysis
#'
#' The Bayesian counterpart of the two-sample (pooled-variance) \emph{t}
#' test. It reports the posterior of the standardized mean difference
#' \eqn{\delta = (\mu_1 - \mu_2)/\sigma} under the default
#' Jeffreys-Zellner-Siow (JZS) prior, a Cauchy prior on \eqn{\delta} with
#' Jeffreys priors on the nuisance parameters (the variance), summarized by
#' its median, mean, a credible interval, and the probability statement
#' \eqn{P(\delta > 0 \mid \mathrm{data})}. The JZS default Bayes factor
#' (Rouder, Speckman, Sun, Morey, & Iverson, 2009) is also reported. See
#' \code{\link{bayes_one_sample_t}} for the model, the package's
#' interpretive stance, and the computational details (exact quadrature,
#' no Monte Carlo error).
#'
#' @param x,y Numeric vectors: the observations of the two independent
#'   groups (\eqn{\delta} is positive when \code{x} runs higher). Omit
#'   both to supply summary statistics instead.
#' @param mean_1,sd_1,n_1 Summary statistics of the first group: mean,
#'   standard deviation, and sample size. The Bayes factor depends on the
#'   data only through the \emph{t} statistic and the sample sizes, so the
#'   summary form is exact, not an approximation. Supply either raw data
#'   or all six summary values, never both.
#' @param mean_2,sd_2,n_2 Summary statistics of the second group.
#' @param prior_location Location of the Cauchy prior on \eqn{\delta}.
#'   Defaults to 0, the JZS prior; a nonzero value centers the prior on an
#'   expected effect (Gronau, Ly, & Wagenmakers, 2020).
#' @param prior_mean,prior_sd Mean and standard deviation of a normal
#'   prior on \eqn{\delta}, for prior beliefs stated as moments. Supplying
#'   them selects the normal prior; they cannot be combined with the
#'   Cauchy arguments. See the prior section of Details.
#' @param prior_scale The way to adjust the prior. It is the scale (width)
#'   of the Cauchy prior on the standardized effect \eqn{\delta} (the JZS
#'   prior), so a user who wants a more or less informative prior sets
#'   \code{prior_scale}: larger values say larger effects are plausible a
#'   priori, smaller values concentrate the prior near zero. The default
#'   \eqn{\sqrt{2}/2 \approx 0.707} is the JZS \dQuote{medium} prior. Fully
#'   custom or subjective priors beyond the Cauchy family are not supported
#'   by the \pkg{BayesFactor} engine.
#' @param conf_level Probability mass of the credible interval.
#'
#' @details
#' The likelihood of the pooled-variance \emph{t} statistic given
#' \eqn{\delta} is noncentral \emph{t} with
#' \eqn{\mathit{df} = n_1 + n_2 - 2} and noncentrality
#' \eqn{\delta \sqrt{n_1 n_2 / (n_1 + n_2)}}; equal variances are assumed,
#' as in the standard JZS development. The raw-scale rows transform the
#' \eqn{\delta} summaries through the pooled standard deviation.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with the same
#'   rows as \code{\link{bayes_one_sample_t}}, plus \code{n_1} and
#'   \code{n_2}.
#'
#'
#' \strong{Specifying the prior.} The default prior on the standardized
#' effect \eqn{\delta} is the JZS Cauchy centered at zero. Its
#' \code{prior_scale} \eqn{r} is not a standard deviation: a Cauchy has
#' no mean and no variance (those integrals diverge), so beliefs stated
#' as prior moments cannot be expressed through it. What the scale does
#' fix is the quartiles: half the prior mass lies within
#' \eqn{\pm r} of the location, so the default \eqn{r = \sqrt{2}/2}
#' says a 50 percent prior bet that \eqn{|\delta| < 0.71}. A directional
#' prior keeps the Cauchy and moves \code{prior_location} (Gronau, Ly, &
#' Wagenmakers, 2020). A researcher who thinks in prior moments instead
#' sets \code{prior_mean} and \code{prior_sd}, which use a normal prior
#' with exactly those moments; the two families are exclusive.
#'
#' The families are also linked by an exact identity: a Cauchy with
#' location \eqn{\mu} and scale \eqn{r} is a normal prior
#' \eqn{N(\mu, r^2/z^2)} whose \eqn{z} is standard normal, that is, a
#' normal prior whose variance you are not sure of. Choosing the Cauchy
#' is therefore choosing a normal prior with built-in doubt about its
#' own width, which is why its tails are heavier and its Bayes factors
#' more conservative. A normal matched to the Cauchy's interquartile
#' range has \code{prior_sd = 1.4826 * prior_scale}. The full posterior
#' of \eqn{\delta} is returned in the \code{"posterior"} attribute as a
#' data frame of \code{delta} and \code{density}, so any posterior
#' probability, not only the reported ones, can be computed from it.
#'
#' A standardized effect size enters through the summary form directly:
#' an observed Cohen's \emph{d} with group sizes \eqn{n_1} and \eqn{n_2}
#' is \code{mean_1 = d, mean_2 = 0, sd_1 = 1, sd_2 = 1}, since \emph{d}
#' is the mean difference in pooled standard deviation units.
#'
#'
#' @references
#' Gronau, Q. F., Ly, A., & Wagenmakers, E.-J. (2020). Informed
#'   Bayesian t-tests. \emph{The American Statistician, 74}(2),
#'   137--143. \doi{10.1080/00031305.2018.1562983}
#'
#' Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., & Iverson, G.
#'   (2009). Bayesian t tests for accepting and rejecting the null
#'   hypothesis. \emph{Psychonomic Bulletin & Review, 16}(2), 225--237.
#'   \doi{10.3758/PBR.16.2.225}
#'
#' Jeffreys, H. (1961). \emph{Theory of probability} (3rd ed.). Oxford
#'   University Press.
#'
#' Zellner, A., & Siow, A. (1980). Posterior odds ratios for selected
#'   regression hypotheses. In J. M. Bernardo, M. H. DeGroot, D. V.
#'   Lindley, & A. F. M. Smith (Eds.), \emph{Bayesian statistics:
#'   Proceedings of the First International Meeting} (pp. 585--603).
#'   University of Valencia Press.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{bayes_one_sample_t}} and
#'   \code{\link{bayes_paired_t}} for the other designs;
#'   \code{\link{ci_smd}} and \code{\link{welch_t}} for the frequentist
#'   analyses of the same comparison.
#'
#' @family Bayesian t analyses
#'
#' @keywords htest
#'
#' @examples
#' set.seed(113)
#' g1 <- rnorm(35, 105, 15)
#' g2 <- rnorm(35, 100, 15)
#' bayes_independent_t(g1, g2)
#'
#' # The probability that the effect is positive is read directly off the
#' # p_delta_positive row: a probability statement, not a p-value.
#'
#' @export
#' @importFrom stats var
bayes_independent_t <- function(x = NULL, y = NULL,
                                mean_1 = NULL, sd_1 = NULL, n_1 = NULL,
                                mean_2 = NULL, sd_2 = NULL, n_2 = NULL,
                                prior_location = 0,
                                prior_scale = sqrt(2) / 2,
                                prior_mean = NULL, prior_sd = NULL,
                                conf_level = 0.95) {
  prior <- .bayes_t_prior(prior_location, prior_scale, prior_mean, prior_sd,
                          cauchy_args_supplied = !missing(prior_scale))
  summaries <- !is.null(mean_1) || !is.null(sd_1) || !is.null(n_1) ||
    !is.null(mean_2) || !is.null(sd_2) || !is.null(n_2)
  if ((is.null(x) || is.null(y)) && !summaries) {
    stop("Supply either raw data through 'x' and 'y' or the summary ",
         "statistics 'mean_1', 'sd_1', 'n_1', 'mean_2', 'sd_2', and ",
         "'n_2'.", call. = FALSE)
  }
  if (!is.null(x) && summaries) {
    stop("Supply raw data through 'x' and 'y' or summary statistics, ",
         "not both.", call. = FALSE)
  }
  if (summaries) {
    .bayes_t_check_summaries(mean = mean_1, sd = sd_1, n = n_1,
                             labels = c("mean_1", "sd_1", "n_1"))
    .bayes_t_check_summaries(mean = mean_2, sd = sd_2, n = n_2,
                             labels = c("mean_2", "sd_2", "n_2"))
  } else {
    .bayes_t_check_x(x, "x")
    .bayes_t_check_x(y, "y")
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }
  if (summaries) {
    n1 <- n_1; n2 <- n_2
    m1 <- mean_1; m2 <- mean_2
    v1 <- sd_1^2; v2 <- sd_2^2
  } else {
    n1 <- length(x); n2 <- length(y)
    m1 <- mean(x); m2 <- mean(y)
    v1 <- stats::var(x); v2 <- stats::var(y)
  }
  df <- n1 + n2 - 2
  sp <- sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / df)
  n_eff <- n1 * n2 / (n1 + n2)
  tt <- (m1 - m2) / (sp * sqrt(1 / n1 + 1 / n2))
  .bayes_t_table(t_obs = tt, df = df, n_eff = n_eff, s_raw = sp,
                 prior = prior, conf_level = conf_level,
                 extra_terms = c("n_1", "n_2"), extra_values = c(n1, n2))
}
