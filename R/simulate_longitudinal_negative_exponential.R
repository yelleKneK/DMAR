# Simulate longitudinal data from the negative exponential growth curve.
#' Simulate Data From a Negative Exponential (Asymptotic Regression) Change Model
#'
#' Generates longitudinal data from a random-coefficients negative
#' exponential change model, also called asymptotic regression (Stevens,
#' 1951): each unit (a person, an animal, a tree) approaches an asymptote at a rate set by a
#' curvature parameter, the parameters vary randomly across units, and
#' each measurement adds level-one error. The deterministic part of the
#' curve is
#' \deqn{\mu(t) = \alpha + \zeta \exp(-\gamma t),}
#' the parameterization of Kelley (2005, 2008). The negative exponential
#' is the simplest of the package's nonlinear change curves: it has one
#' asymptote and no point of inflection, so it describes change that is
#' fastest at the first assessment and decelerates thereafter.
#'
#' @param n A single positive integer, the number of units (persons,
#'   animals, trees, classrooms) whose trajectories are drawn, or a
#'   vector giving the number of units for each population, one entry
#'   per parameter vector in \code{fixed_parameters}.
#' @param target_times Numeric vector of the nominal measurement times,
#'   one shared schedule for every unit. Give either this or
#'   \code{time_range}.
#' @param time_range Alternative to \code{target_times}: \code{c(lower,
#'   upper)} bounds from which each unit draws its own measurement
#'   times, so no two units share a schedule (e.g., age in weeks at
#'   testing rather than a fixed grade). Requires \code{occasions};
#'   with \code{time_range}, the level-one error must be a single
#'   \code{error_variance} with the default independent structure, and
#'   \code{timing_sd} does not apply.
#' @param occasions With \code{time_range}: a single positive integer
#'   (every unit measured the same number of times) or \code{c(min,
#'   max)}, from which each unit's number of measurement times is drawn
#'   uniformly.
#' @param time_distribution Distribution of the unit-specific times over
#'   \code{time_range}; currently \code{"uniform"}.
#' @param fixed_parameters The population parameters
#'   \code{c(alpha = , zeta = , gamma = )} (an unnamed length-3 vector is
#'   taken in that order), or a list of such
#'   vectors, one per population (distinct parameter vectors):
#'   \describe{
#'     \item{\code{alpha}}{The asymptote: the value \eqn{\mu(t)}
#'       approaches as \eqn{t} grows.}
#'     \item{\code{zeta}}{The negative of the total change: the curve
#'       starts at the intercept \eqn{\phi = \alpha + \zeta} and travels
#'       \eqn{-\zeta} units to the asymptote. Positive \code{zeta} gives
#'       asymptotic decay toward \code{alpha} from above; negative
#'       \code{zeta} gives asymptotic growth from below.}
#'     \item{\code{gamma}}{The curvature (\eqn{\gamma > 0}): the rate at
#'       which the remaining distance to the asymptote closes per unit
#'       time. Larger values reach the asymptote sooner.}
#'   }
#' @param random_variances Between-unit variances of
#'   \code{(alpha, zeta, gamma)}, a single number recycled to all three
#'   or a length-3 vector. Default \code{0} (a common curve for
#'   everyone). A named vector over any subset of the parameter
#'   names (e.g., \code{c(gamma = 0.02)}) varies only those named and
#'   leaves the rest fixed.
#' @param random_correlation Optional 3-by-3 correlation matrix among the
#'   random parameters; default uncorrelated.
#' @param error_variance,reliability,error_structure,error_correlation,timing_sd
#'   The level-one error and assessment-time machinery, with the same
#'   meaning as in \code{\link{simulate_longitudinal_polynomial}}:
#'   specify exactly one of \code{error_variance} (scalar, per-occasion
#'   vector, or full covariance matrix) or \code{reliability} (a target
#'   average per-occasion reliability, solved by the delta method here);
#'   \code{error_structure} and \code{error_correlation} set the
#'   across-occasion error correlation; \code{timing_sd} jitters each
#'   unit's actual assessment times around the nominal targets.
#'
#' @return A long-format \code{data.frame} with columns \code{id},
#'   \code{population}, \code{occasion}, \code{target_time}, \code{time},
#'   \code{true_score}, and \code{y}, directly usable with
#'   \code{\link{plot_trajectories}} and nonlinear mixed-model fitters
#'   such as \code{nlme::nlme()}. Attributes carry the \code{model}, the
#'   \code{fixed_parameters}, the between-unit covariance
#'   \code{random_covariance}, the level-one \code{error_variance} and
#'   \code{error_covariance}, and \code{reliability_by_occasion} (from
#'   the first-order delta method true-score variance, which can drift
#'   from the realized variance ratio when the random variances are
#'   large relative to the mean curve). The \code{schedule} attribute
#'   records \code{"shared"} or \code{"unit_specific"}.
#'
#' @details
#' Every parameter answers a substantive question: where does change end
#' (\code{alpha}), where does it start (\eqn{\phi = \alpha + \zeta}),
#' and how fast does the gap close (\code{gamma})? That
#' interpretability is the argument for nonlinear change models over
#' polynomials, whose coefficients describe no landmark of the process
#' (Kelley, 2005, 2008); the package vignette on nonlinear growth
#' develops the comparison.
#'
#' @references
#' Kelley, K. (2005). \emph{Estimating nonlinear change models in
#'   heterogeneous populations when class membership is unknown: Defining
#'   and developing the latent classification differential change model}
#'   (Doctoral dissertation). University of Notre Dame.
#'
#' Kelley, K. (2008). Nonlinear change models in populations with
#'   unobserved heterogeneity. \emph{Methodology, 4}(3), 97--112.
#'
#' Stevens, W. L. (1951). Asymptotic regression. \emph{Biometrics,
#'   7}(3), 247--267.
#'
#' @seealso \code{\link{simulate_longitudinal_logistic}},
#'   \code{\link{simulate_longitudinal_gompertz}},
#'   \code{\link{simulate_longitudinal_richards}} for the sigmoidal
#'   members of the family; \code{\link{simulate_longitudinal_polynomial}}
#'   for the polynomial counterpart; \code{\link{plot_trajectories}} for
#'   plotting the result.
#'
#' @examples
#' # The six-curve illustration from Kelley (2005): three growth curves
#' # (intercept 0, asymptote 1) that differ only in curvature, and
#' # three decay curves (intercept 1, asymptote 0) that mirror them.
#' # Each curve is its own population of size one, so a single call draws the
#' # whole panel.
#' panel <- simulate_longitudinal_negative_exponential(
#'   n = 1, target_times = seq(0, 10, by = 0.1),
#'   fixed_parameters = list(
#'     c(alpha = 1, zeta = -1, gamma = 0.9),
#'     c(alpha = 1, zeta = -1, gamma = 0.4),
#'     c(alpha = 1, zeta = -1, gamma = 0.2),
#'     c(alpha = 0, zeta =  1, gamma = 1.2),
#'     c(alpha = 0, zeta =  1, gamma = 0.5),
#'     c(alpha = 0, zeta =  1, gamma = 0.3)),
#'   error_variance = 0
#' )
#' plot_trajectories(panel, id = "id", time = "time",
#'                   outcome = "true_score", group = "population")
#'
#' # Individual differences in a single parameter: only the curvature
#' # varies (a named entry leaves the other variances at zero), so all
#' # trajectories share their start and their destination but close the
#' # gap at their own rates.
#' set.seed(113)
#' d_gamma <- simulate_longitudinal_negative_exponential(
#'   n = 25, target_times = seq(0, 8, by = 0.5),
#'   fixed_parameters = c(alpha = 100, zeta = -80, gamma = 0.5),
#'   random_variances = c(gamma = 0.02), error_variance = 0
#' )
#' plot_trajectories(d_gamma, id = "id", time = "time",
#'                   outcome = "true_score")
#'
#' # Individual differences in every parameter, plus level-one error:
#' # vocabulary learning that starts near 20 words (phi = alpha + zeta),
#' # climbs toward an asymptote near 100, and closes about 40% of the
#' # remaining gap per month (gamma = 0.5).
#' set.seed(113)
#' d <- simulate_longitudinal_negative_exponential(
#'   n = 30, target_times = 0:8,
#'   fixed_parameters = c(alpha = 100, zeta = -80, gamma = 0.5),
#'   random_variances = c(alpha = 25, zeta = 16, gamma = 0.01),
#'   error_variance = 9
#' )
#' head(d)
#' plot_trajectories(d, id = "id", time = "time", outcome = "y")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords datagen
#'
#' @family data simulators
#'
#' @export
simulate_longitudinal_negative_exponential <- function(n,
                                                       target_times = NULL,
                                                       fixed_parameters,
                                                       time_range = NULL,
                                                       occasions = NULL,
                                                       time_distribution = "uniform",
                                           random_variances = 0,
                                                       random_correlation = NULL,
                                                       error_variance = NULL,
                                                       reliability = NULL,
                                                       error_structure = c("independent",
                                                                           "ar1",
                                                                           "compound_symmetry",
                                                                           "toeplitz"),
                                                       error_correlation = NULL,
                                                       timing_sd = 0) {
  .simulate_longitudinal_growth(
    model = "negative_exponential",
    n = n, target_times = target_times,
    fixed_parameters = fixed_parameters,
    time_range = time_range, occasions = occasions,
    time_distribution = time_distribution,
    random_variances = random_variances,
    random_correlation = random_correlation,
    error_variance = error_variance, reliability = reliability,
    error_structure = match.arg(error_structure),
    error_correlation = error_correlation, timing_sd = timing_sd
  )
}
