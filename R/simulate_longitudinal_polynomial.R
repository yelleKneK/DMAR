#' Simulate Data From a Polynomial Change (Growth) Model
#'
#' Generates longitudinal data from a random-coefficients polynomial change
#' model: each subject follows a degree-\eqn{P} polynomial in time whose
#' coefficients vary randomly across subjects, and each measurement adds
#' independent level-one error. The polynomial order is general (order 0 is a
#' flat line, 1 a straight line, 2 a quadratic, and so on), one or several
#' populations may differ in their mean trajectories, the level-one error can be set
#' directly or pinned to a target measurement reliability, and the actual time
#' of each assessment may jitter away from its nominal target. This is the
#' Monte Carlo companion to \code{\link{ss_power_pcm}}, which plans power for
#' the same model under the closed-form assumptions of Raudenbush and Liu
#' (2001); the simulator can relax those assumptions (notably the assumption of
#' fixed, error-free, equally reliable assessment times) and study what happens.
#'
#' @param n A single positive integer (equal number of units in every population) or a
#'   numeric vector of length \code{G} giving the number of subjects in each of
#'   the \code{G} populations, where \code{G} is the number of mean trajectories
#'   supplied through \code{fixed_coefficients}.
#' @param target_times A numeric vector of the nominal (planned) measurement
#'   times, length \eqn{M}, one shared schedule for every unit. They need
#'   not be equally spaced. A degree-\eqn{P} model requires
#'   \eqn{M \ge P + 1} occasions. Give either this or \code{time_range}.
#' @param time_range Alternative to \code{target_times}: \code{c(lower,
#'   upper)} bounds from which each unit draws its own measurement times,
#'   so no two units share a schedule (e.g., age in weeks at testing
#'   rather than a fixed grade). Requires \code{occasions}; with
#'   \code{time_range}, the level-one error must be a single
#'   \code{error_variance} with the default independent structure, and
#'   \code{timing_sd} does not apply.
#' @param occasions With \code{time_range}: a single positive integer
#'   (every unit measured the same number of times) or \code{c(min,
#'   max)}, from which each unit's number of measurement times is drawn
#'   uniformly. Every value must be at least \eqn{P + 1} so each unit's
#'   trajectory identifies the polynomial.
#' @param time_distribution Distribution of the unit-specific times over
#'   \code{time_range}; currently \code{"uniform"}.
#' @param fixed_coefficients The population mean trajectory, as the coefficients
#'   of a polynomial in time, ordered from the intercept upward:
#'   \code{c(b0, b1, ..., bP)} encodes
#'   \eqn{b_0 + b_1 t + b_2 t^2 + \dots + b_P t^P}. The length sets the
#'   polynomial order \eqn{P} (length 1 is order 0, a flat line at \code{b0}).
#'   For several populations, pass a \emph{list} of equal-length coefficient
#'   vectors, one per population; the populations then differ in their mean
#'   trajectories but share
#'   the variance components below (the Raudenbush-Liu two-group setup).
#' @param random_variances The between-subject variances of the polynomial
#'   coefficients (the diagonal of the level-two covariance matrix), as a single
#'   number recycled to all \eqn{P + 1} coefficients or a vector of length
#'   \eqn{P + 1}. Default \code{0}, a fixed common trajectory with no
#'   between-subject heterogeneity. Entries may be zero coefficient by
#'   coefficient (e.g., a random intercept with a fixed slope).
#' @param random_correlation Optional \eqn{(P + 1) \times (P + 1)} correlation
#'   matrix among the random coefficients. Default \code{NULL} treats them as
#'   uncorrelated. Combined with \code{random_variances} it forms the level-two
#'   covariance \eqn{T = D R D}, with \eqn{D} the diagonal matrix of
#'   coefficient standard deviations.
#' @param error_variance The level-one (within-subject) measurement error
#'   variance \eqn{\sigma^2_e}. Most simply a single non-negative number (the
#'   same error variance at every occasion). It may instead be a length-\eqn{M}
#'   vector of per-occasion (heteroscedastic) error variances, or a full
#'   \eqn{M \times M} error covariance matrix when the errors are correlated
#'   across occasions in an arbitrary (unstructured) way. For the common
#'   structured cases, give a scalar or vector here and set
#'   \code{error_structure} and \code{error_correlation} instead of building the
#'   matrix by hand. Specify exactly one of \code{error_variance} or
#'   \code{reliability}.
#' @param reliability A target measurement reliability in \eqn{(0, 1)} from
#'   which \eqn{\sigma^2_e} is derived (see Details). The value is the
#'   \emph{average} per-occasion reliability across \code{target_times}; the
#'   per-occasion reliabilities, which generally differ from one another, are
#'   returned in the \code{"reliability_by_occasion"} attribute. Requires at
#'   least one positive entry in \code{random_variances}. The solved error
#'   variance is homoscedastic and may still be given an across-occasion
#'   correlation through \code{error_structure}. Specify exactly one of
#'   \code{error_variance} or \code{reliability}.
#' @param error_structure The correlation pattern of the level-one errors across
#'   occasions: \code{"independent"} (the default, uncorrelated errors),
#'   \code{"ar1"} (a first-order autoregressive decay \eqn{\rho^{|j-k|}}, so
#'   errors at occasions closer in time are more alike), \code{"compound_symmetry"}
#'   (a constant correlation \eqn{\rho} between every pair of occasions), or
#'   \code{"toeplitz"} (a banded structure set by the lag-1 through
#'   lag-\eqn{(M-1)} correlations). Ignored when \code{error_variance} is a full
#'   covariance matrix, which already fixes the structure.
#' @param error_correlation The correlation parameter(s) for
#'   \code{error_structure}: a single number \eqn{\rho} for \code{"ar1"} and
#'   \code{"compound_symmetry"}, or a vector of the lag-1 to lag-\eqn{(M-1)}
#'   correlations for \code{"toeplitz"}. Left \code{NULL} for
#'   \code{"independent"}. The implied correlation matrix must be positive
#'   semidefinite.
#' @param timing_sd The standard deviation of the difference between a subject's
#'   actual and nominal assessment time, as a single number recycled to all
#'   occasions or a vector of length \eqn{M}. Default \code{0}, every subject
#'   measured exactly on schedule. A positive value draws each subject's actual
#'   time at occasion \eqn{m} as \eqn{\tau_m + N(0, \text{timing\_sd}_m^2)} and
#'   evaluates that subject's true score at the actual time, while the nominal
#'   target is retained in a separate column (see Details).
#'
#' @return A long-format \code{data.frame} with one row per subject-occasion and
#'   the columns
#'   \describe{
#'     \item{\code{id}}{Factor uniquely identifying each subject.}
#'     \item{\code{population}}{Factor with \code{G} levels (\code{"1"}, ...) giving
#'       each unit's population (its data generating parameter vector). With one parameter vector there is one level.}
#'     \item{\code{occasion}}{Integer occasion index, \code{1} to \eqn{M}.}
#'     \item{\code{target_time}}{The nominal (planned) measurement time.}
#'     \item{\code{time}}{The actual measurement time (equal to
#'       \code{target_time} when \code{timing_sd = 0}, otherwise jittered).}
#'     \item{\code{true_score}}{The subject's latent trajectory value at the
#'       actual time, before level-one error.}
#'     \item{\code{y}}{The observed score, \code{true_score} plus level-one
#'       error.}
#'   }
#'   The returned object carries attributes \code{"error_variance"} (the
#'   \eqn{\sigma^2_e} used, a scalar when the errors are homoscedastic and
#'   independent, otherwise the vector of per-occasion error variances),
#'   \code{"error_covariance"} (the full \eqn{M \times M} level-one error
#'   covariance actually used), \code{"reliability_by_occasion"} (the
#'   per-occasion reliabilities at the nominal times),
#'   \code{"random_covariance"} (the level-two covariance \eqn{T}),
#'   \code{"polynomial_order"} (\eqn{P}), and \code{"schedule"}
#'   (\code{"shared"} or \code{"unit_specific"}). With \code{time_range}
#'   there is no shared occasion grid, so \code{"error_covariance"} is
#'   \code{NA} and \code{"reliability_by_occasion"} is \code{NA}. The
#'   format is directly usable with
#'   \code{\link{plot_trajectories}} and with mixed-model fitters such as
#'   \code{nlme::lme()} or \code{lme4::lmer()}.
#'
#' @details
#' \strong{The model.} Subject \eqn{i} in population \eqn{g} has a random coefficient
#' vector \eqn{\pi_i = (\pi_{i0}, \dots, \pi_{iP})} drawn from a multivariate
#' normal with mean the population's \code{fixed_coefficients} \eqn{\beta_g} and
#' covariance \eqn{T}. The latent trajectory is the polynomial
#' \eqn{\mu_i(t) = \sum_{k=0}^{P} \pi_{ik}\, t^k}, and the observed score at a
#' measurement time \eqn{t} is \eqn{y = \mu_i(t) + e}, with
#' \eqn{e \sim N(0, \sigma^2_e)} independent across occasions. Order 0 collapses
#' to a flat line \eqn{\mu_i(t) = \pi_{i0}}; order 1 is the straight-line growth
#' model underlying Raudenbush and Liu (2001) and \code{\link{ss_power_pcm}}.
#'
#' \strong{Coefficient metric.} The coefficients here are the ordinary (raw)
#' polynomial coefficients on \eqn{t^k}, which is the most transparent metric
#' for specifying a trajectory. The derivative-scaled change coefficient used by
#' \code{\link{ss_power_pcm}} and the Raudenbush-Liu power formulas is \eqn{P!}
#' times the leading (highest-order) coefficient supplied here, so a quadratic
#' with \code{fixed_coefficients = c(b0, b1, b2)} corresponds to a
#' Raudenbush-Liu quadratic change coefficient of \eqn{2! \, b_2 = 2 b_2}.
#'
#' \strong{Reliability varies by occasion.} At a measurement time \eqn{t} the
#' implied between-subject (true-score) variance is the quadratic form
#' \eqn{c(t)^\top T\, c(t)} with \eqn{c(t) = (1, t, t^2, \dots, t^P)^\top}, so
#' the classical reliability of the observed score,
#' \deqn{\rho_{XX}(t) = \frac{c(t)^\top T\, c(t)}{c(t)^\top T\, c(t) + \sigma^2_e},}
#' generally \emph{changes from occasion to occasion}: a growth measurement is
#' not equally reliable everywhere, because the spread of true scores depends on
#' where in time you measure relative to the centering of the polynomial and the
#' random-effect structure. When \code{reliability} is supplied, \eqn{\sigma^2_e}
#' is solved (by \code{\link[stats]{uniroot}}) so that the average of
#' \eqn{\rho_{XX}(t)} over the nominal \code{target_times} equals the requested
#' value; the occasion-by-occasion reliabilities are returned in the
#' \code{"reliability_by_occasion"} attribute so the variation is visible rather
#' than hidden behind a single number. Reliability is only meaningful when there
#' is true-score variance to detect, so this route requires
#' \code{random_variances} to be positive for at least one coefficient.
#'
#' \strong{Measurement errors need not be independent or equal.} The simplest
#' model adds an independent, equal-variance error at every occasion, but
#' repeated measurements of the same person are often correlated (an unmodeled
#' state, a rater, or an instrument carries over from one wave to the next) and
#' may be more or less variable at different waves. The level-one errors are
#' drawn from \eqn{N(0, \Sigma_e)}, and \eqn{\Sigma_e} can be set three ways: a
#' scalar or per-occasion \code{error_variance} combined with an
#' \code{error_structure} (\code{"ar1"} for autoregressive decay, the natural
#' choice when occasions are ordered in time; \code{"compound_symmetry"} for an
#' equicorrelated error; \code{"toeplitz"} for a general banded pattern), or a
#' full covariance matrix passed directly as \code{error_variance}. Because
#' classical reliability at an occasion is a marginal quantity, it depends only
#' on the diagonal of \eqn{\Sigma_e}; the across-occasion error correlation
#' leaves \code{"reliability_by_occasion"} unchanged but does affect how a
#' mixed model that assumes independent errors performs, which is exactly the
#' kind of misspecification this simulator is meant to let a user study.
#'
#' \strong{Assessment timing is rarely exact.} Designs are written as if every
#' subject is measured at the same fixed times (\dQuote{the 7-day follow-up}),
#' but in practice people arrive early or late, so the actual time differs from
#' the nominal target. Setting \code{timing_sd > 0} draws each subject's actual
#' time per occasion and evaluates the true score \emph{at the time the
#' measurement really happened}, while \code{target_time} keeps the nominal
#' value an analyst would typically use. Analyzing on the nominal time when the
#' data were in fact collected on jittered times biases estimates of the change
#' coefficients, and the bias grows with the order of the trend and with the
#' size of the timing variability. The two time columns let a user quantify that
#' bias by fitting the same model on \code{time} versus \code{target_time}.
#'
#' \strong{Why the closed-form Raudenbush-Liu planner does not cover all of
#' this.} The power formulas in Raudenbush and Liu (2001), carried by
#' \code{\link{ss_power_pcm}}, are exact under three assumptions this simulator
#' can relax: every subject is measured at the \emph{same}, equally spaced,
#' \emph{error-free} occasion times; the level-one error variance is a single
#' constant (so the closed form needs no notion of an occasion-varying
#' reliability); and the within-subject sampling variance of the change
#' coefficient has the known form \eqn{V = \sigma^2_e f^{2p} (M - p - 1)! /
#' [K_p (M + p)!]}. Those assumptions buy a clean formula, but real designs
#' violate them: assessments drift in time, and reliability is not the same at
#' every wave. This function is the Monte Carlo complement that lets a
#' researcher generate data under the messier reality and check how far the
#' closed-form power and the fitted estimates can be trusted.
#'
#' @references
#' Kelley, K., & Rausch, J. R. (2011). Sample size planning for longitudinal
#'   models: Accuracy in parameter estimation for polynomial change parameters.
#'   \emph{Psychological Methods, 16}(4), 391--405. \doi{10.1037/a0023352}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 15 on the analysis of repeated measures
#'   and growth.)
#'
#' Raudenbush, S. W., & Liu, X.-F. (2001). Effects of study duration, frequency
#'   of observation, and sample size on power in studies of group differences in
#'   polynomial change. \emph{Psychological Methods, 6}(4), 387--401.
#'   \doi{10.1037/1082-989X.6.4.387}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_power_pcm}} for closed-form power planning on the
#'   same model, \code{\link{plot_trajectories}} to visualize the simulated
#'   curves, and \code{\link[MASS]{mvrnorm}} for the random-coefficient draw.
#'
#' @family data simulators
#'
#' @keywords design datagen
#'
#' @examples
#' # 1. One population of linear growers with a random intercept and a random slope,
#' #    measured yearly for four years (five occasions), with the level-one
#' #    error set directly.
#' set.seed(113)
#' d <- simulate_longitudinal_polynomial(
#'   n                  = 50,
#'   target_times       = 0:4,
#'   fixed_coefficients = c(10, 1.5),          # intercept 10, slope 1.5 per year
#'   random_variances   = c(4, 0.25),          # var(intercept) = 4, var(slope) = .25
#'   error_variance     = 1
#' )
#' head(d)
#'
#' # 2. Two populations that differ only in their slope (a treatment that changes
#' #    the rate of growth). Pass a list of coefficient vectors, one per population.
#' set.seed(113)
#' two <- simulate_longitudinal_polynomial(
#'   n                  = c(40, 40),
#'   target_times       = 0:4,
#'   fixed_coefficients = list(control = c(10, 1.0), treatment = c(10, 1.8)),
#'   random_variances   = c(4, 0.25),
#'   error_variance     = 1
#' )
#' aggregate(y ~ population + occasion, data = two, FUN = mean)
#'
#' # 3. Pin the level-one error to a target reliability instead of setting it
#' #    directly. The single number is the average reliability across occasions;
#' #    the per-occasion values differ and are returned as an attribute.
#' set.seed(113)
#' rel <- simulate_longitudinal_polynomial(
#'   n                  = 100,
#'   target_times       = 0:4,
#'   fixed_coefficients = c(10, 1.5),
#'   random_variances   = c(4, 0.25),
#'   reliability        = 0.80
#' )
#' attr(rel, "error_variance")
#' round(attr(rel, "reliability_by_occasion"), 3)   # not constant across waves
#'
#' # 4. Assessment-time jitter: the nominal "yearly" schedule, but subjects
#' #    actually arrive a little early or late (SD of about six weeks on a
#' #    one-year scale). The nominal and actual times are kept in separate
#' #    columns so the consequences of analyzing on the nominal time can be
#' #    studied.
#' set.seed(113)
#' jit <- simulate_longitudinal_polynomial(
#'   n                  = 30,
#'   target_times       = 0:4,
#'   fixed_coefficients = c(10, 1.5),
#'   random_variances   = c(4, 0.25),
#'   error_variance     = 1,
#'   timing_sd          = 0.12
#' )
#' head(jit[, c("id", "occasion", "target_time", "time")])
#'
#' # 4b. Unit-specific measurement times: each of 12 children is tested
#' #     between 40 and 90 weeks of age, five to nine times, no two on the
#' #     same schedule. The level-one error is a single variance; the
#' #     "schedule" attribute records the design.
#' set.seed(113)
#' ages <- simulate_longitudinal_polynomial(
#'   n                  = 12,
#'   time_range         = c(40, 90),
#'   occasions          = c(5, 9),
#'   fixed_coefficients = c(10, 0.5),
#'   random_variances   = c(4, 0.01),
#'   error_variance     = 2
#' )
#' attr(ages, "schedule")
#' table(table(ages$id))   # units per occasion count
#'
#' # 5. A flat line (order 0): no growth, only a random subject level and
#' #    measurement error. The coefficient vector has length one.
#' set.seed(113)
#' flat <- simulate_longitudinal_polynomial(
#'   n                  = 20,
#'   target_times       = 0:4,
#'   fixed_coefficients = 5,
#'   random_variances   = 2,
#'   error_variance     = 1
#' )
#' head(flat)
#'
#' # 6. Autocorrelated measurement error: the same error variance at each wave,
#' #    but the level-one errors decay as an AR(1) process (errors at adjacent
#' #    occasions correlate 0.5), the kind of dependence a model assuming
#' #    independent errors would miss. The full error covariance is returned.
#' set.seed(113)
#' ar <- simulate_longitudinal_polynomial(
#'   n                  = 40,
#'   target_times       = 0:4,
#'   fixed_coefficients = c(10, 1.5),
#'   random_variances   = c(4, 0.25),
#'   error_variance     = 1,
#'   error_structure    = "ar1",
#'   error_correlation  = 0.5
#' )
#' round(attr(ar, "error_covariance"), 3)
#'
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm uniroot
simulate_longitudinal_polynomial <- function(n,
                                              target_times = NULL,
                                              fixed_coefficients,
                                              time_range = NULL,
                                              occasions = NULL,
                                              time_distribution = "uniform",
                                              random_variances   = 0,
                                              random_correlation = NULL,
                                              error_variance     = NULL,
                                              reliability        = NULL,
                                              error_structure    = c("independent",
                                                                      "ar1",
                                                                      "compound_symmetry",
                                                                      "toeplitz"),
                                              error_correlation  = NULL,
                                              timing_sd          = 0) {

  error_structure <- match.arg(error_structure)

  # ---------- Resolve populations and polynomial order from fixed_coefficients ----------
  if (is.list(fixed_coefficients)) {
    if (length(fixed_coefficients) < 1L) {
      stop("'fixed_coefficients' must contain at least one coefficient vector.",
           call. = FALSE)
    }
    ok_elt <- vapply(fixed_coefficients,
                     function(b) is.numeric(b) && length(b) >= 1L && !anyNA(b),
                     logical(1))
    if (!all(ok_elt)) {
      stop("Each element of 'fixed_coefficients' must be a non-empty numeric ",
           "vector with no missing values.", call. = FALSE)
    }
    lens <- vapply(fixed_coefficients, length, integer(1))
    if (any(lens != lens[1L])) {
      stop("Every population's coefficient vector in 'fixed_coefficients' must have ",
           "the same length (the same polynomial order across populations).",
           call. = FALSE)
    }
    coef_list <- lapply(fixed_coefficients, as.numeric)
  } else {
    if (!is.numeric(fixed_coefficients) || length(fixed_coefficients) < 1L ||
        anyNA(fixed_coefficients)) {
      stop("'fixed_coefficients' must be a numeric vector (one population) or a ",
           "list of equal-length numeric vectors (several populations), with no ",
           "missing values.", call. = FALSE)
    }
    coef_list <- list(as.numeric(fixed_coefficients))
  }
  G <- length(coef_list)
  P <- length(coef_list[[1L]]) - 1L     # polynomial order (0 = flat line)

  # ---------- Measurement schedule ----------
  if (is.null(target_times) == is.null(time_range)) {
    stop("Specify exactly one of 'target_times' (one shared schedule for ",
         "every unit) or 'time_range' (unit-specific times drawn between ",
         "a lower and an upper bound).", call. = FALSE)
  }
  individual_times <- !is.null(time_range)
  if (individual_times) {
    time_distribution <- match.arg(time_distribution, c("uniform"))
    if (!is.numeric(time_range) || length(time_range) != 2L ||
        any(!is.finite(time_range)) || time_range[1L] >= time_range[2L]) {
      stop("'time_range' must be c(lower, upper) with lower < upper.",
           call. = FALSE)
    }
    if (is.null(occasions) || !is.numeric(occasions) || anyNA(occasions) ||
        !length(occasions) %in% c(1L, 2L) || any(occasions < 1) ||
        any(occasions != round(occasions)) ||
        (length(occasions) == 2L && occasions[1L] > occasions[2L])) {
      stop("With 'time_range', give 'occasions' as a single positive ",
           "integer (the same number of measurement times for every unit) ",
           "or c(min, max) for a unit-varying number.", call. = FALSE)
    }
    occasions <- as.integer(occasions)
    if (min(occasions) < P + 1L) {
      stop(sprintf(paste0("A degree-%d polynomial needs at least %d ",
                          "measurement occasions per unit, but 'occasions' ",
                          "allows as few as %d. Raise 'occasions' or lower ",
                          "the polynomial order set by 'fixed_coefficients'."),
                   P, P + 1L, min(occasions)), call. = FALSE)
    }
    M <- max(occasions)
  } else {
    if (!is.numeric(target_times) || length(target_times) < 1L ||
        any(!is.finite(target_times))) {
      stop("'target_times' must be a numeric vector of finite measurement times.",
           call. = FALSE)
    }
    target_times <- as.numeric(target_times)
    M <- length(target_times)
    if (M < P + 1L) {
      stop(sprintf(paste0("A degree-%d polynomial needs at least %d measurement ",
                          "occasions, but 'target_times' has %d. Add occasions or ",
                          "lower the polynomial order set by 'fixed_coefficients'."),
                   P, P + 1L, M), call. = FALSE)
    }
  }

  # ---------- Per-population unit counts ----------
  if (length(n) == 1L) {
    if (!is.numeric(n) || is.na(n) || n < 1 || n != round(n)) {
      stop("'n' must be a single positive integer (units per population).",
           call. = FALSE)
    }
    n_per_population <- rep(as.integer(n), G)
  } else if (length(n) == G) {
    if (!is.numeric(n) || anyNA(n) || any(n < 1) || any(n != round(n))) {
      stop("Each entry of 'n' must be a positive integer.", call. = FALSE)
    }
    n_per_population <- as.integer(n)
  } else {
    stop(sprintf(paste0("'n' must be a single positive integer or a numeric ",
                        "vector of length %d (one count of units per population)."),
                 G), call. = FALSE)
  }

  # ---------- Level-two (between-subject) covariance ----------
  if (!is.numeric(random_variances) || anyNA(random_variances) ||
      any(random_variances < 0)) {
    stop("'random_variances' must be non-negative.", call. = FALSE)
  }
  if (length(random_variances) == 1L) {
    random_variances <- rep(as.numeric(random_variances), P + 1L)
  } else if (length(random_variances) != P + 1L) {
    stop(sprintf(paste0("'random_variances' must have length 1 or %d (one ",
                        "variance per polynomial coefficient, orders 0 to %d)."),
                 P + 1L, P), call. = FALSE)
  }

  if (is.null(random_correlation)) {
    Rmat <- diag(P + 1L)
  } else {
    Rmat <- as.matrix(random_correlation)
    if (nrow(Rmat) != P + 1L || ncol(Rmat) != P + 1L) {
      stop(sprintf("'random_correlation' must be a %d-by-%d matrix.",
                   P + 1L, P + 1L), call. = FALSE)
    }
    if (any(abs(Rmat - t(Rmat)) > 1e-8)) {
      stop("'random_correlation' must be symmetric.", call. = FALSE)
    }
    if (any(abs(diag(Rmat) - 1) > 1e-8)) {
      stop("'random_correlation' must have 1 on the diagonal (it is a ",
           "correlation matrix).", call. = FALSE)
    }
    ev <- eigen(Rmat, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < -1e-8) {
      stop("'random_correlation' must be positive semidefinite.", call. = FALSE)
    }
  }
  std_dev <- sqrt(random_variances)
  Tcov    <- outer(std_dev, std_dev) * Rmat     # diag(std_dev) %*% Rmat %*% diag(std_dev)

  # ---------- Assessment-time jitter ----------
  if (!is.numeric(timing_sd) || anyNA(timing_sd) || any(timing_sd < 0)) {
    stop("'timing_sd' must be non-negative.", call. = FALSE)
  }
  if (individual_times) {
    if (!(length(timing_sd) == 1L && timing_sd == 0)) {
      stop("'timing_sd' does not combine with 'time_range': unit-specific ",
           "times are already drawn per unit, so jitter around a nominal ",
           "schedule has no meaning there.", call. = FALSE)
    }
  } else if (length(timing_sd) == 1L) {
    timing_sd <- rep(as.numeric(timing_sd), M)
  } else if (length(timing_sd) != M) {
    stop(sprintf(paste0("'timing_sd' must have length 1 or %d (one standard ",
                        "deviation per occasion)."), M), call. = FALSE)
  }

  # ---------- Level-one error covariance: from error_variance or reliability ----------
  # The within-subject errors are drawn from N(0, Sigma_e), an M-by-M
  # covariance. Sigma_e comes from one of three routes: (a) error_variance as a
  # scalar (homoscedastic) or length-M vector (heteroscedastic) of marginal
  # variances combined with the error_structure correlation pattern; (b)
  # error_variance as a full M-by-M covariance matrix (an unstructured Sigma);
  # or (c) reliability, which solves a single homoscedastic error variance so
  # the average per-occasion reliability hits the target and then applies the
  # chosen correlation structure on top. Per-occasion (classical) reliability
  # depends only on the marginal error variances (the diagonal of Sigma_e), so
  # the across-occasion correlation never changes the reported reliabilities.
  if (is.null(error_variance) == is.null(reliability)) {
    stop("Specify exactly one of 'error_variance' or 'reliability'.",
         call. = FALSE)
  }
  if (individual_times) {
    # With unit-specific times there is no shared occasion grid, so the
    # per-occasion machinery (reliability targets, heteroscedastic or
    # correlated errors) has nothing to attach to.
    if (!is.null(reliability)) {
      stop("'reliability' needs a shared measurement schedule: ",
           "per-occasion reliability is undefined when every unit has its ",
           "own times. Give 'error_variance' instead.", call. = FALSE)
    }
    if (!is.numeric(error_variance) || length(error_variance) != 1L ||
        is.na(error_variance) || error_variance < 0) {
      stop("With 'time_range', 'error_variance' must be a single ",
           "non-negative number: per-occasion vectors and covariance ",
           "matrices need a shared schedule.", call. = FALSE)
    }
    if (error_structure != "independent" || !is.null(error_correlation)) {
      stop("With 'time_range', 'error_structure' must be \"independent\": ",
           "across-occasion error correlation needs a shared schedule.",
           call. = FALSE)
    }
    err_scalar <- as.numeric(error_variance)
    Sigma_e <- NULL
    err_var_by_occ <- err_scalar
    reliability_by_occ <- NA_real_
    err_indep <- TRUE
  } else {
  # Implied true-score variance c(t)'T c(t) at each nominal occasion.
  basis_target     <- outer(target_times, 0:P, `^`)            # M x (P + 1)
  true_var_by_occ  <- rowSums((basis_target %*% Tcov) * basis_target)

  if (!is.null(error_variance) && is.matrix(error_variance)) {
    # (b) Full covariance matrix: error_variance IS Sigma_e.
    Sigma_e <- error_variance
    if (nrow(Sigma_e) != M || ncol(Sigma_e) != M) {
      stop(sprintf(paste0("When 'error_variance' is a covariance matrix it must ",
                          "be %d-by-%d (one row and column per occasion)."),
                   M, M), call. = FALSE)
    }
    if (!is.numeric(Sigma_e) || anyNA(Sigma_e)) {
      stop("The 'error_variance' covariance matrix must be numeric with no ",
           "missing values.", call. = FALSE)
    }
    if (any(abs(Sigma_e - t(Sigma_e)) > 1e-8)) {
      stop("The 'error_variance' covariance matrix must be symmetric.",
           call. = FALSE)
    }
    if (min(eigen(Sigma_e, symmetric = TRUE, only.values = TRUE)$values) < -1e-8) {
      stop("The 'error_variance' covariance matrix must be positive ",
           "semidefinite.", call. = FALSE)
    }
    if (error_structure != "independent" || !is.null(error_correlation)) {
      stop("When 'error_variance' is a full covariance matrix it already ",
           "specifies the error structure; leave 'error_structure' and ",
           "'error_correlation' at their defaults.", call. = FALSE)
    }
    err_var_by_occ <- diag(Sigma_e)
  } else {
    if (!is.null(error_variance)) {
      # (a) Scalar or per-occasion variances, plus a correlation structure.
      if (!is.numeric(error_variance) || anyNA(error_variance) ||
          any(error_variance < 0) ||
          !(length(error_variance) %in% c(1L, M))) {
        stop(sprintf(paste0("'error_variance' must be a single non-negative ",
                            "number, a length-%d vector of per-occasion ",
                            "variances, or a %d-by-%d covariance matrix."),
                     M, M, M), call. = FALSE)
      }
      err_var_by_occ <- if (length(error_variance) == 1L)
        rep(as.numeric(error_variance), M) else as.numeric(error_variance)
    } else {
      # (c) Solve a single homoscedastic error variance from a target average
      # reliability, then apply the correlation structure below.
      if (!is.numeric(reliability) || length(reliability) != 1L ||
          is.na(reliability) || reliability <= 0 || reliability >= 1) {
        stop("'reliability' must be a single number in (0, 1).", call. = FALSE)
      }
      if (!any(random_variances > 0)) {
        stop("'reliability' requires positive between-subject variance: set ",
             "'random_variances' above 0 for at least one coefficient, or ",
             "specify 'error_variance' directly.", call. = FALSE)
      }
      # As sigma2_e -> 0+, each occasion with positive true-score variance has
      # reliability -> 1 and each with zero true-score variance stays at 0, so
      # the largest achievable average reliability is the fraction of occasions
      # with positive true-score variance.
      pos      <- true_var_by_occ > 0
      sup_rel  <- mean(pos)
      if (reliability >= sup_rel) {
        stop(sprintf(paste0("The requested average reliability (%.3f) is not ",
                            "achievable for this design: %d of %d nominal ",
                            "occasions have positive true-score variance, so ",
                            "the average per-occasion reliability cannot reach ",
                            "%.3f even with zero measurement error. Lower ",
                            "'reliability', or change the design or ",
                            "'random_variances'."),
                     reliability, sum(pos), M, sup_rel), call. = FALSE)
      }
      rel_gap <- function(s) sum(true_var_by_occ[pos] /
                                   (true_var_by_occ[pos] + s)) / M - reliability
      upper <- max(true_var_by_occ[pos])
      while (rel_gap(upper) > 0) upper <- upper * 2
      sigma2_e <- stats::uniroot(rel_gap, lower = 0, upper = upper,
                                 tol = .Machine$double.eps^0.5)$root
      err_var_by_occ <- rep(sigma2_e, M)
    }
    R_err   <- .error_cor_matrix(error_structure, error_correlation, M)
    sd_occ  <- sqrt(err_var_by_occ)
    Sigma_e <- outer(sd_occ, sd_occ) * R_err
  }

  # Errors are independent across occasions exactly when Sigma_e is diagonal.
  err_indep <- max(abs(Sigma_e - diag(diag(Sigma_e), M))) < 1e-12
  denom <- true_var_by_occ + err_var_by_occ
  reliability_by_occ <- true_var_by_occ / denom
  reliability_by_occ[denom == 0] <- NA_real_   # 0/0 when both variances are 0
  names(reliability_by_occ) <- paste0("occasion_", seq_len(M))
  }

  # ---------- Generate ----------
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required. Install with install.packages(\"MASS\").",
         call. = FALSE)
  }
  any_random <- any(random_variances > 0)
  N_total    <- sum(n_per_population)

  if (individual_times) {
    # Each unit draws its own number of occasions and its own sorted times;
    # rows accumulate per unit rather than into a preallocated grid.
    rows_all <- vector("list", G)
    subj <- 0L
    for (g in seq_len(G)) {
      beta_g <- coef_list[[g]]
      n_g    <- n_per_population[g]
      coefs_g <- if (any_random) {
        matrix(MASS::mvrnorm(n_g, mu = beta_g, Sigma = Tcov), nrow = n_g)
      } else {
        matrix(beta_g, nrow = n_g, ncol = P + 1L, byrow = TRUE)
      }
      rows <- vector("list", n_g)
      for (i in seq_len(n_g)) {
        subj <- subj + 1L
        m_i <- if (length(occasions) == 2L) {
          counts <- seq.int(occasions[1L], occasions[2L])
          counts[sample.int(length(counts), 1L)]
        } else {
          occasions
        }
        t_i     <- sort(stats::runif(m_i, time_range[1L], time_range[2L]))
        basis_i <- outer(t_i, 0:P, `^`)
        true_i  <- as.numeric(basis_i %*% coefs_g[i, ])
        y_i <- if (err_scalar > 0) {
          true_i + stats::rnorm(m_i, 0, sqrt(err_scalar))
        } else {
          true_i
        }
        rows[[i]] <- data.frame(
          id          = subj,
          population  = g,
          occasion    = seq_len(m_i),
          target_time = t_i,
          time        = t_i,
          true_score  = true_i,
          y           = y_i,
          stringsAsFactors = FALSE
        )
      }
      rows_all[[g]] <- do.call(rbind, rows)
    }
    tmp <- do.call(rbind, rows_all)
    out <- data.frame(
      id          = factor(tmp$id, levels = seq_len(N_total)),
      population  = factor(tmp$population, levels = seq_len(G),
                           labels = as.character(seq_len(G))),
      occasion    = tmp$occasion,
      target_time = tmp$target_time,
      time        = tmp$time,
      true_score  = tmp$true_score,
      y           = tmp$y,
      stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
  } else {

  any_jitter <- any(timing_sd > 0)
  sd_occ_e   <- sqrt(diag(Sigma_e))    # per-occasion error SD (the diagonal)
  draw_error <- any(diag(Sigma_e) > 0)
  total_rows <- N_total * M

  id_vec     <- integer(total_rows)
  population_vec  <- integer(total_rows)
  occ_vec    <- integer(total_rows)
  target_vec <- numeric(total_rows)
  time_vec   <- numeric(total_rows)
  true_vec   <- numeric(total_rows)
  y_vec      <- numeric(total_rows)

  row  <- 1L
  subj <- 0L
  for (g in seq_len(G)) {
    beta_g <- coef_list[[g]]
    n_g    <- n_per_population[g]

    if (any_random) {
      coefs_g <- matrix(MASS::mvrnorm(n_g, mu = beta_g, Sigma = Tcov), nrow = n_g)
    } else {
      coefs_g <- matrix(beta_g, nrow = n_g, ncol = P + 1L, byrow = TRUE)
    }

    # Correlated errors are drawn for the whole population at once; independent
    # errors (the default) stay on the per-subject rnorm() path below so the
    # homoscedastic case reproduces the simpler model's random stream exactly.
    if (draw_error && !err_indep) {
      E_g <- matrix(MASS::mvrnorm(n_g, mu = rep(0, M), Sigma = Sigma_e),
                    nrow = n_g)
    }

    for (i in seq_len(n_g)) {
      subj <- subj + 1L
      t_actual <- if (any_jitter) target_times + stats::rnorm(M, 0, timing_sd) else target_times
      basis_i  <- outer(t_actual, 0:P, `^`)            # M x (P + 1)
      true_i   <- as.numeric(basis_i %*% coefs_g[i, ])
      if (!draw_error) {
        y_i <- true_i
      } else if (err_indep) {
        y_i <- true_i + stats::rnorm(M, 0, sd_occ_e)
      } else {
        y_i <- true_i + E_g[i, ]
      }

      idx              <- row:(row + M - 1L)
      id_vec[idx]      <- subj
      population_vec[idx]   <- g
      occ_vec[idx]     <- seq_len(M)
      target_vec[idx]  <- target_times
      time_vec[idx]    <- t_actual
      true_vec[idx]    <- true_i
      y_vec[idx]       <- y_i
      row <- row + M
    }
  }

  out <- data.frame(
    id          = factor(id_vec, levels = seq_len(N_total)),
    population  = factor(population_vec, levels = seq_len(G),
                         labels = as.character(seq_len(G))),
    occasion    = occ_vec,
    target_time = target_vec,
    time        = time_vec,
    true_score  = true_vec,
    y           = y_vec,
    stringsAsFactors = FALSE
  )
  }

  # "error_variance" stays a scalar in the common homoscedastic, independent
  # case (backward compatible); otherwise it is the vector of per-occasion
  # variances. The full M-by-M error covariance is available in
  # "error_covariance" whenever the schedule is shared.
  attr(out, "error_variance") <- if (individual_times) {
    err_scalar
  } else if (err_indep && max(abs(err_var_by_occ - err_var_by_occ[1L])) < 1e-12) {
    err_var_by_occ[1L]
  } else {
    err_var_by_occ
  }
  attr(out, "error_covariance")        <- if (individual_times) NA else Sigma_e
  attr(out, "reliability_by_occasion") <- reliability_by_occ
  attr(out, "random_covariance")       <- Tcov
  attr(out, "polynomial_order")        <- P
  attr(out, "schedule") <- if (individual_times) "unit_specific" else "shared"
  out
}


# Build an M-by-M level-one error correlation matrix for a named structure,
# validating the supplied correlation parameter(s). Not exported; the one
# consumer is simulate_longitudinal_polynomial(). "independent" is the identity;
# "ar1" decays as rho^|j - k| (measurements closer in time more alike);
# "compound_symmetry" is a constant off-diagonal rho; "toeplitz" takes a vector
# of lag-1 to lag-(M - 1) correlations. The result is checked to be positive
# semidefinite so the draw is from a legitimate covariance.
.error_cor_matrix <- function(structure, rho, M) {
  if (structure == "independent") {
    if (!is.null(rho))
      stop("'error_correlation' is only used when 'error_structure' is not ",
           "\"independent\".", call. = FALSE)
    return(diag(M))
  }
  if (is.null(rho))
    stop(sprintf("'error_structure = \"%s\"' requires 'error_correlation'.",
                 structure), call. = FALSE)
  if (!is.numeric(rho) || anyNA(rho))
    stop("'error_correlation' must be numeric with no missing values.",
         call. = FALSE)
  lag <- abs(outer(seq_len(M), seq_len(M), "-"))
  R <- switch(structure,
    ar1 = {
      if (length(rho) != 1L)
        stop("AR(1) 'error_correlation' must be a single number.", call. = FALSE)
      if (abs(rho) >= 1)
        stop("AR(1) 'error_correlation' must lie in (-1, 1).", call. = FALSE)
      rho ^ lag
    },
    compound_symmetry = {
      if (length(rho) != 1L)
        stop("compound symmetry 'error_correlation' must be a single number.",
             call. = FALSE)
      if (abs(rho) > 1)
        stop("compound symmetry 'error_correlation' must lie in [-1, 1].",
             call. = FALSE)
      R0 <- matrix(rho, M, M); diag(R0) <- 1; R0
    },
    toeplitz = {
      if (length(rho) != M - 1L)
        stop(sprintf(paste0("Toeplitz 'error_correlation' must have length %d ",
                            "(the lag-1 to lag-%d correlations)."),
                     M - 1L, M - 1L), call. = FALSE)
      if (any(abs(rho) > 1))
        stop("Toeplitz 'error_correlation' values must lie in [-1, 1].",
             call. = FALSE)
      matrix(c(1, rho)[lag + 1L], M, M)
    }
  )
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < -1e-8)
    stop(sprintf(paste0("The implied '%s' error correlation matrix is not ",
                        "positive semidefinite; choose a smaller ",
                        "'error_correlation'."), structure), call. = FALSE)
  R
}
